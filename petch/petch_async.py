import datetime
from colorama import Fore, Style
import gzip
import os
import time
import re
import pandas as pd
import aiohttp
import asyncio
from io import StringIO
from tqdm import tqdm
import argparse

# Async function to fetch and decompress data
async def fetch_data(session, url):
    async with session.get(url) as response:
        response.raise_for_status()
        compressed_data = await response.read()
        return gzip.decompress(compressed_data).decode('utf-8')

# Main async function to fetch Uniprot proteome entries
async def fetch_uniprot_proteomes_async(filterout=None):
    ref_url = "https://rest.uniprot.org/proteomes/stream?compressed=true&fields=upid%2Corganism%2Cprotein_count%2Clineage&format=tsv&query=%28*%29+AND+%28proteome_type%3A1%29"
    other_url = "https://rest.uniprot.org/proteomes/stream?compressed=true&fields=upid%2Corganism%2Cprotein_count%2Clineage&format=tsv&query=%28*%29+AND+%28proteome_type%3A2%29&sort=cpd+asc"

    try:
        async with aiohttp.ClientSession() as session:
            print("Downloading Reference Proteomes...")
            ref_task = fetch_data(session, ref_url)
            print("Downloading Other Proteomes...")
            other_task = fetch_data(session, other_url)

            # Run both tasks concurrently
            ref_data, other_data = await asyncio.gather(ref_task, other_task)

            # Process Reference Proteomes
            ref_proteomes_df = pd.read_csv(StringIO(ref_data), sep='\t').dropna(subset=['Organism'])
            print(f"{ref_proteomes_df.shape[0]} Reference Proteome Entries.")

            if filterout:
                for keyword in filterout:
                    ref_proteomes_df = ref_proteomes_df[~ref_proteomes_df['Taxonomic lineage'].fillna('').str.contains(keyword)]
                    print(f"{ref_proteomes_df.shape[0]} Reference Entries after filtering for '{keyword}'.")

            # Process Other Proteomes
            other_proteomes_df = pd.read_csv(StringIO(other_data), sep='\t').dropna(subset=['Organism'])
            print(f"{other_proteomes_df.shape[0]} Other Proteome Entries.")

            if filterout:
                for keyword in filterout:
                    other_proteomes_df = other_proteomes_df[~other_proteomes_df['Taxonomic lineage'].fillna('').str.contains(keyword)]
                    print(f"{other_proteomes_df.shape[0]} Other Entries after filtering for '{keyword}'.")

            return ref_proteomes_df, other_proteomes_df
            
    except aiohttp.ClientError as e:
        print(f"Network-related error occurred: {e}")
    except asyncio.CancelledError:
        print("The task was canceled.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Wrapper function to run the async code
def fetch_uniprot_proteomes(filterout=None):
    return asyncio.run(fetch_uniprot_proteomes_async(filterout=filterout))

# Read the taxa list from a CSV or TSV file
def read_taxa(input_file):
    try:
        # Check file extension to determine delimiter
        file_extension = os.path.splitext(input_file)[1].lower()
        if file_extension == '.csv':
            df = pd.read_csv(input_file)
        elif file_extension == '.tsv':
            df = pd.read_csv(input_file, sep='\t')
        else:
            raise ValueError("Unsupported file format. Please provide a CSV or TSV file.")
        
        # Display information about the loaded file
        original_rows, num_columns = df.shape
        print(f"{input_file} loaded: {original_rows} Rows, {num_columns} Column(s)")

        # Process taxa column
        taxa_names = df['taxa'].str.lower().drop_duplicates().tolist()
        
        # Check for and handle duplicates
        removed_duplicates = original_rows - len(taxa_names)
        if removed_duplicates > 0:
            print(f"{removed_duplicates} duplicate(s) removed.")
        else:
            print(f"No duplicates found.")
        
        return taxa_names

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        return None
    except Exception as e:
        print(f"Error reading file '{input_file}': {e}")
        return None

# Find best match from the ref_proteomes or other_proteomes
def find_all_matches(taxa_list, ref_proteomes, other_proteomes, output_dir):
    matched_df = pd.DataFrame(columns=ref_proteomes.columns)
    unmatched_taxa = []

    for taxa in taxa_list:
        escaped_taxa = re.escape(taxa)
        ref_match = ref_proteomes[ref_proteomes['Organism'].str.contains(escaped_taxa, case=False)]
        
        if ref_match.empty:
            other_match = other_proteomes[other_proteomes['Organism'].str.contains(escaped_taxa, case=False)]
            if not other_match.empty:
                other_match = other_match.copy()
                other_match.loc[:, 'Proteome Source'] = 'Other'
                other_match.loc[:, 'Queried Term'] = taxa
                matched_df = pd.concat([matched_df, other_match.drop(columns=['Taxonomic lineage']).iloc[[0]]])
            else:
                unmatched_taxa.append(taxa)
        else:
            ref_match = ref_match.copy()
            ref_match.loc[:, 'Proteome Source'] = 'Ref'
            ref_match.loc[:, 'Queried Term'] = taxa
            matched_df = pd.concat([matched_df, ref_match.drop(columns=['Taxonomic lineage']).iloc[[0]]])

    matched_df = matched_df.drop_duplicates(subset=['Proteome Id'], keep='first')

    os.makedirs(output_dir, exist_ok=True)
    if unmatched_taxa:
        unmatched_df = pd.DataFrame({'Unmatched Taxa': unmatched_taxa})
        unmatched_df.to_csv(os.path.join(output_dir, "no_matches.csv"), index=False)
        print(f"Unmatched taxa saved to '{output_dir}/no_matches.csv'.")

    matched_df = matched_df.dropna(axis=1, how='all')
    matched_df.to_csv(os.path.join(output_dir, "matches.csv"), index=False)
    print(f"Matched data saved to '{output_dir}/matches.csv'")

    return matched_df


# Download proteomes of matched taxa from Uniprot
async def download_proteome_async(proteome_id, taxa):
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=proteome:{proteome_id}"
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as response:
                if response.status == 200:
                    data = await response.text()
                    return data, taxa
                else:
                    return None, taxa
    except aiohttp.ClientError as e:
        print(f"Network-related error occurred: {e}")
        return None, taxa
    except Exception as e:
        print(f"Unexpected error occurred: {e}")
        return None, taxa

def fetch_fastas(matched_df, output_dir):
    total_rows = len(matched_df)
    output_file = os.path.join(output_dir, "Petch.fasta")

    # Remove existing FASTA file if it exists
    if os.path.exists(output_file):
        os.remove(output_file)
        print(f"Existing file '{output_file}' has been removed.")

    print(f"Fetching {total_rows} proteomes...")

    async def fetch_fasta_data():
        tasks = []
        progress_bar = tqdm(total=total_rows, desc="Downloading Proteomes", unit="proteome")
        
        # Loop through the DataFrame and create tasks for each row
        for index, row in matched_df.iterrows():
            proteome_id = row['Proteome Id']
            taxa = row['Queried Term']
            task = asyncio.create_task(download_proteome_async(proteome_id, taxa))
            task.add_done_callback(lambda _: progress_bar.update())
            tasks.append(task)

        results = await asyncio.gather(*tasks)
        progress_bar.close()
        return results

    # Run the async function and process results
    results = asyncio.run(fetch_fasta_data())

    # Write to the fasta file based on the results
    with open(output_file, "a") as fasta_file:
        for proteome_data, taxa in results:
            if proteome_data:
                fasta_file.write(proteome_data)
            else:
                print(f"Error downloading proteome for '{taxa}'")




#Count the number of entries per organism in FASTA.
def count_entries_per_organism(fasta_file, output_dir):
    organism_counts = {}
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'): 
                organism = line.split('OS=')[1].split(' OX=')[0]
                organism_counts[organism] = organism_counts.get(organism, 0) + 1

    organism_counts_df = pd.DataFrame(organism_counts.items(), columns=['Fasta Entry', 'Entry Count'])
    organism_counts_df.to_csv(os.path.join(output_dir, "fasta_counts.csv"), index=False)
    return organism_counts_df



#Validation using Fasta Entires and Proteome Counts.
def merge_dataframes(matched_df, organism_counts_df, output_dir):
    preprocess_organism = lambda text: ' '.join(text.split()[:2]).rstrip('.')
    matched_df['Organism'] = matched_df['Organism'].apply(preprocess_organism)
    organism_counts_df['Fasta Entry'] = organism_counts_df['Fasta Entry'].apply(preprocess_organism)
   
    merged_df = pd.merge(matched_df, organism_counts_df, left_on='Organism', right_on='Fasta Entry', how='right')
    merged_df.fillna({'Fasta Entry Count': 'N/A'}, inplace=True)
    merged_df = merged_df[['Organism', 'Fasta Entry', 'Protein count', 'Entry Count']]
    merged_df['Counts Match'] = merged_df.apply(lambda row: 'pass' if row['Protein count'] == row['Entry Count'] else '-', axis=1)
    merged_df.to_csv(os.path.join(output_dir, "valid_counts.csv"), index=False)
    
    print(f"Validation Merge completed, saved to '{output_dir}/valid_counts.csv'.")

    if merged_df['Counts Match'].eq('pass').all():
        print("Validation - PASSED -")
    else:
        print("Validation - FAILED -")

    return merged_df


def parse_arguments():
    parser = argparse.ArgumentParser(description="Fetch UniProt proteomes and filter by taxa list")
    parser.add_argument("input_file", help="CSV or TSV file containing taxa list")
    parser.add_argument("--filterout", nargs='+', help="Keywords to filter out")
    parser.add_argument("--output_dir", default="outputs", help="Directory to save output files (default: 'outputs')")
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    start_time = time.time()
    args = parse_arguments()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    print(f"{Fore.YELLOW}** Running PETCH: Proteome Fetcher Script **{Style.RESET_ALL}")

    print(f"{Fore.BLUE}Hold On A Sec, Grapping Data From Uniprot...{Style.RESET_ALL}")
    ref_proteomes, other_proteomes = fetch_uniprot_proteomes(args.filterout)

    print(f"{Fore.BLUE}Inspecting Taxa...{Style.RESET_ALL}")
    taxa_list = read_taxa(args.input_file)

    print(f"{Fore.BLUE}Looking For Matches...{Style.RESET_ALL}")
    matched_df = find_all_matches(taxa_list, ref_proteomes, other_proteomes, output_dir)

    print(f"{Fore.BLUE}Fetching FASTAs...{Style.RESET_ALL}")
    fetch_fastas(matched_df, output_dir)

    print(f"{Fore.BLUE}Counting Entries Per Organism...{Style.RESET_ALL}")
    organism_counts_df = count_entries_per_organism(os.path.join(output_dir, "Petch.fasta"), output_dir)

    print(f"{Fore.BLUE}Comparing Dataframes...{Style.RESET_ALL}")
    merge_dataframes(matched_df, organism_counts_df, output_dir)

    print(f"{Fore.GREEN}Script Execution Completed, Proteomes Acquired!!{Style.RESET_ALL}")
    print(f"Results saved to: {output_dir}")