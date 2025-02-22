import datetime
import gzip
import os
import time
import re
import pandas as pd
import aiohttp
import asyncio
import math
from io import StringIO
from tqdm import tqdm
import argparse
from colorama import Fore, Style
import concurrent.futures

# Constants
TIMEOUT = aiohttp.ClientTimeout(total=60)  # 60 seconds
REF_URL = "https://rest.uniprot.org/proteomes/stream?compressed=true&fields=upid%2Corganism%2Cprotein_count%2Clineage&format=tsv&query=%28*%29+AND+%28proteome_type%3A1%29"
OTHER_URL = "https://rest.uniprot.org/proteomes/stream?compressed=true&fields=upid%2Corganism%2Cprotein_count%2Clineage&format=tsv&query=%28*%29+AND+%28proteome_type%3A2%29&sort=cpd+asc"

async def fetch_data(session, url, retries=2):
    """
    Fetch and decompress data from a given URL asynchronously.

    Args:
        session (aiohttp.ClientSession): The aiohttp session.
        url (str): The URL to fetch data from.
        retries (int): Number of retry attempts.

    Returns:
        str: Decompressed data as a string.
    """
    for attempt in range(retries):
        try:
            async with session.get(url) as response:
                response.raise_for_status()
                compressed_data = await response.read()
                return gzip.decompress(compressed_data).decode('utf-8')
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            if attempt == retries - 1:
                raise
            await asyncio.sleep(2 ** attempt)  # Exponential backoff

async def fetch_uniprot_proteomes_async(filterout=None):
    """
    Fetch UniProt proteomes (reference and other) asynchronously.

    Args:
        filterout (list, optional): Keywords to filter out from the proteomes.

    Returns:
        tuple: Two DataFrames containing reference and other proteomes.
    """
    try:
        async with aiohttp.ClientSession(timeout=TIMEOUT) as session:
            print("Downloading Reference Proteomes...")
            ref_task = fetch_data(session, REF_URL)
            print("Downloading Other Proteomes...")
            other_task = fetch_data(session, OTHER_URL)

            ref_data, other_data = await asyncio.gather(ref_task, other_task)

            ref_proteomes_df = pd.read_csv(StringIO(ref_data), sep='\t').dropna(subset=['Organism'])
            print(f"{ref_proteomes_df.shape[0]} Reference Proteome Entries.")

            if filterout:
                ref_proteomes_df = filter_proteomes(ref_proteomes_df, filterout, "Reference")

            other_proteomes_df = pd.read_csv(StringIO(other_data), sep='\t').dropna(subset=['Organism'])
            print(f"{other_proteomes_df.shape[0]} Other Proteome Entries.")

            if filterout:
                other_proteomes_df = filter_proteomes(other_proteomes_df, filterout, "Other")

            return ref_proteomes_df, other_proteomes_df
            
    except aiohttp.ClientError as e:
        print(f"Network-related error occurred: {e}")
    except asyncio.CancelledError:
        print("The task was canceled.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def filter_proteomes(df, keywords, proteome_type):
    """
    Filter proteomes DataFrame based on keywords.

    Args:
        df (pd.DataFrame): The DataFrame to filter.
        keywords (list): Keywords to filter out.
        proteome_type (str): Type of proteome (e.g., "Reference" or "Other").

    Returns:
        pd.DataFrame: Filtered DataFrame.
    """
    for keyword in keywords:
        df = df[~df['Taxonomic lineage'].fillna('').str.contains(keyword)]
        print(f"{df.shape[0]} {proteome_type} Entries after filtering for '{keyword}'.")
    return df

def fetch_uniprot_proteomes(filterout=None):
    """
    Wrapper function to run the async proteome fetching code.

    Args:
        filterout (list, optional): Keywords to filter out.

    Returns:
        tuple: Two DataFrames containing reference and other proteomes.
    """
    return asyncio.run(fetch_uniprot_proteomes_async(filterout=filterout))

def read_taxa(input_file):
    """
    Read taxa list from a CSV or TSV file.

    Args:
        input_file (str): Path to the input file.

    Returns:
        list: List of unique taxa names.
    """
    try:
        file_extension = os.path.splitext(input_file)[1].lower()
        delimiter = ',' if file_extension == '.csv' else '\t' if file_extension == '.tsv' else None
        if delimiter is None:
            raise ValueError("Unsupported file format. Please provide a CSV or TSV file.")
        
        df = pd.read_csv(input_file, sep=delimiter)
        original_rows, num_columns = df.shape
        print(f"{input_file} loaded: {original_rows} Rows, {num_columns} Column(s)")

        taxa_names = df['taxa'].str.lower().drop_duplicates().tolist()
        removed_duplicates = original_rows - len(taxa_names)
        print(f"{removed_duplicates} duplicate(s) removed." if removed_duplicates > 0 else "No duplicates found.")
        
        return taxa_names

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        return None
    except Exception as e:
        print(f"Error reading file '{input_file}': {e}")
        return None

def find_all_matches(taxa_list, ref_proteomes, other_proteomes, output_dir):
    """
    Find matches for taxa in reference and other proteomes.

    Args:
        taxa_list (list): List of taxa to match.
        ref_proteomes (pd.DataFrame): Reference proteomes DataFrame.
        other_proteomes (pd.DataFrame): Other proteomes DataFrame.
        output_dir (str): Directory to save output files.

    Returns:
        pd.DataFrame: DataFrame containing matched proteomes.
    """
    matched_df = pd.DataFrame(columns=ref_proteomes.columns)
    unmatched_taxa = []

    def match_taxa(taxa):
        escaped_taxa = re.escape(taxa)
        ref_match = ref_proteomes[ref_proteomes['Organism'].str.contains(escaped_taxa, case=False)]
        
        if ref_match.empty:
            other_match = other_proteomes[other_proteomes['Organism'].str.contains(escaped_taxa, case=False)]
            if not other_match.empty:
                other_match = other_match.copy()
                other_match.loc[:, 'Proteome Source'] = 'Other'
                other_match.loc[:, 'Queried Term'] = taxa
                return other_match.drop(columns=['Taxonomic lineage']).iloc[[0]], None
            else:
                return None, taxa
        else:
            ref_match = ref_match.copy()
            ref_match.loc[:, 'Proteome Source'] = 'Ref'
            ref_match.loc[:, 'Queried Term'] = taxa
            return ref_match.drop(columns=['Taxonomic lineage']).iloc[[0]], None

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(match_taxa, taxa) for taxa in taxa_list]
        for future in concurrent.futures.as_completed(futures):
            match, unmatched = future.result()
            if match is not None:
                matched_df = pd.concat([matched_df, match])
            if unmatched is not None:
                unmatched_taxa.append(unmatched)

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

async def download_proteome_async(session, proteome_id, taxa, retries=3):
    """
    Download a proteome from UniProt asynchronously.

    Args:
        session (aiohttp.ClientSession): The aiohttp session.
        proteome_id (str): The proteome ID to download.
        taxa (str): The taxa name for logging.
        retries (int): Number of retry attempts.

    Returns:
        tuple: (Decompressed data, taxa) or (None, taxa) if failed.
    """
    url = f"https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=proteome:{proteome_id}"
    
    for attempt in range(retries):
        try:
            async with session.get(url) as response:
                if response.status == 200:
                    compressed_data = await response.read() 
                    decompressed_data = gzip.decompress(compressed_data).decode("utf-8")
                    return decompressed_data, taxa
                else:
                    print(f"Attempt {attempt + 1} failed for {taxa}: HTTP {response.status}")
                    if attempt == retries - 1:
                        return None, taxa
                    await asyncio.sleep(2 ** attempt)  # Exponential backoff
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            print(f"Attempt {attempt + 1} failed for {taxa}: {e}")
            if attempt == retries - 1:
                return None, taxa
            await asyncio.sleep(2 ** attempt)  # Exponential backoff

async def fetch_fasta_data(matched_df, output_file, batch_size=50):
    """
    Fetch FASTA data for matched proteomes in batches.

    Args:
        matched_df (pd.DataFrame): DataFrame containing matched proteomes.
        output_file (str): Path to the output FASTA file.
        batch_size (int): Number of proteomes to fetch per batch.
    """
    async with aiohttp.ClientSession(timeout=TIMEOUT) as session:
        total_rows = len(matched_df)
        progress_bar = tqdm(total=total_rows, desc="Downloading Proteomes", unit="proteome")
        
        num_batches = math.ceil(total_rows / batch_size)
        for batch_num in range(num_batches):
            start_idx = batch_num * batch_size
            end_idx = start_idx + batch_size
            batch_df = matched_df.iloc[start_idx:end_idx]

            tasks = []
            for index, row in batch_df.iterrows():
                proteome_id = row['Proteome Id']
                taxa = row['Queried Term']
                task = asyncio.create_task(download_proteome_async(session, proteome_id, taxa))
                task.add_done_callback(lambda _: progress_bar.update())
                tasks.append(task)

            results = await asyncio.gather(*tasks)
            with open(output_file, "a") as fasta_file:
                for proteome_data, taxa in results:
                    if proteome_data:
                        fasta_file.write(proteome_data)
                    else:
                        print(f"Error downloading proteome for '{taxa}'")

            await asyncio.sleep(1)  # Rate limiting

        progress_bar.close()

def fetch_fastas(matched_df, output_dir, batch_size=50):
    """
    Fetch FASTA files for matched proteomes.

    Args:
        matched_df (pd.DataFrame): DataFrame containing matched proteomes.
        output_dir (str): Directory to save the output FASTA file.
        batch_size (int): Number of proteomes to fetch per batch.
    """
    output_file = os.path.join(output_dir, "Petch.fasta")
    if os.path.exists(output_file):
        os.remove(output_file)
        print(f"Existing file '{output_file}' has been removed.")

    print(f"Fetching {len(matched_df)} proteomes in batches of {batch_size}...")
    asyncio.run(fetch_fasta_data(matched_df, output_file, batch_size))

def count_entries_per_organism(fasta_file, output_dir):
    """
    Count the number of FASTA entries per organism.

    Args:
        fasta_file (str): Path to the input FASTA file.
        output_dir (str): Directory to save the output counts file.

    Returns:
        pd.DataFrame: DataFrame containing organism counts.
    """
    organism_counts = {}
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'): 
                organism = line.split('OS=')[1].split(' OX=')[0]
                organism_counts[organism] = organism_counts.get(organism, 0) + 1

    organism_counts_df = pd.DataFrame(organism_counts.items(), columns=['Fasta Entry', 'Entry Count'])
    organism_counts_df.to_csv(os.path.join(output_dir, "fasta_counts.csv"), index=False)
    return organism_counts_df

def merge_dataframes(matched_df, organism_counts_df, output_dir):
    """
    Merge matched proteomes with organism counts for validation.

    Args:
        matched_df (pd.DataFrame): DataFrame containing matched proteomes.
        organism_counts_df (pd.DataFrame): DataFrame containing organism counts.
        output_dir (str): Directory to save the validation results.

    Returns:
        pd.DataFrame: Merged DataFrame with validation results.
    """
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
    return parser.parse_args()

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
    fetch_fastas(matched_df, output_dir, batch_size=50)

    print(f"{Fore.BLUE}Counting Entries Per Organism...{Style.RESET_ALL}")
    organism_counts_df = count_entries_per_organism(os.path.join(output_dir, "Petch.fasta"), output_dir)

    print(f"{Fore.BLUE}Comparing Dataframes...{Style.RESET_ALL}")
    merge_dataframes(matched_df, organism_counts_df, output_dir)

    print(f"{Fore.GREEN}Script Execution Completed, Proteomes Acquired!!{Style.RESET_ALL}")
    print(f"Results saved to: {output_dir}")
    print(f"{Fore.YELLOW}Total Execution Time: {(time.time() - start_time) / 60:.2f} minutes{Style.RESET_ALL}")