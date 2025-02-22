import aiohttp
import asyncio
import pandas as pd
import re
import gzip
import os
import time
import argparse
from io import StringIO
from bs4 import BeautifulSoup
from tqdm import tqdm
from colorama import Fore, Style

# Constants
BASE_URL = "https://www.ncbi.nlm.nih.gov/research/cog/cogcategory/"
COG_API_URL = "https://www.ncbi.nlm.nih.gov/research/cog/api/cogdef/?cog="
FILTERED_TAXA_FILE = 'filtered_uniprot_data.tsv'

def normalize_gene_name(gene_name):
    """Normalize gene names by removing numeric suffixes."""
    if pd.isna(gene_name):
        return None
    return re.sub(r'[_\d]+$', '', str(gene_name))

async def fetch_cog_gene_names(session, cog_id):
    """Fetch Gene Names for a given COG ID from NCBI API."""
    async with session.get(f"{COG_API_URL}{cog_id}", timeout=10) as response:
        response.raise_for_status()
        data = await response.json()
        if "results" in data and data["results"]:
            return [normalize_gene_name(gene) for gene in data["results"][0].get("genes", [])]
    return []

async def fetch_page(session, url):
    """Fetch a single page of results asynchronously and parse HTML."""
    async with session.get(url) as response:
        content = await response.text()
        soup = BeautifulSoup(content, 'html.parser')
        return soup

async def fetch_category_data(session, category):
    """Fetch and scrape data for a single COG category."""
    url = BASE_URL + category
    soup = await fetch_page(session, url)
    table = soup.find('table', {'id': 'cog_result_table'})

    if not table:
        print(f"Table with id 'cog_result_table' not found for category {category}.")
        return None

    headers = [header.get_text(strip=True) for header in table.find_all('th')]
    all_results = [ [col.get_text(strip=True) for col in row.find_all('td')] 
                   for row in table.find_all('tr')[1:] ]

    return pd.DataFrame(all_results, columns=headers)

async def fetch_all_categories(categories):
    """Fetch and scrape all COG categories asynchronously."""
    async with aiohttp.ClientSession() as session:
        all_data = [await fetch_category_data(session, category) for category in categories]
        return pd.concat([df for df in all_data if df is not None], ignore_index=True)

async def fetch_gene_data(session, gene, output_file, retries=3, batch_size=100):
    """Fetch gene data from UniProt's API for a single gene with pagination."""
    query_url = f"https://rest.uniprot.org/uniprotkb/search?compressed=true&fields=accession%2Cid%2Cprotein_name%2Corganism_name%2Cgene_primary&format=tsv&query=gene:{gene}&size=500"
    batch_data = []

    try:
        while query_url:
            async with session.get(query_url) as response:
                response.raise_for_status()
                decompressed_data = gzip.decompress(await response.read()).decode('utf-8')
                df = pd.read_csv(StringIO(decompressed_data), sep='\t')
                batch_data.append(df)

                if len(batch_data) >= batch_size:
                    pd.concat(batch_data, ignore_index=True).to_csv(output_file, mode='a', header=not os.path.exists(output_file), index=False)
                    batch_data = []

                next_page_url = response.headers.get("link")
                query_url = next_page_url.split(";")[0].strip("<>") if next_page_url and 'rel="next"' in next_page_url else None

        if batch_data:
            pd.concat(batch_data, ignore_index=True).to_csv(output_file, mode='a', header=not os.path.exists(output_file), index=False)

        return True

    except Exception as e:
        if retries > 0:
            await asyncio.sleep(2 ** (3 - retries))
            return await fetch_gene_data(session, gene, output_file, retries - 1)
        print(f"Failed to fetch data for gene {gene}: {e}")
        return False

async def fetch_all_genes(cog_genes, output_file, max_retries=3):
    """Fetch data for all COG genes concurrently using asynchronous requests."""
    if os.path.exists(output_file):
        os.remove(output_file)
        print(f"Existing output file {output_file} removed.")

    async with aiohttp.ClientSession() as session:
        failed_genes = []
        retry_count = 0

        while retry_count <= max_retries:
            if retry_count > 0:
                print(f"Retrying failed genes (attempt {retry_count}/{max_retries})...")
                genes_to_fetch = failed_genes.copy()
                failed_genes = []
            else:
                genes_to_fetch = cog_genes

            tasks = [fetch_gene_data(session, gene, output_file) for gene in genes_to_fetch]
            print("Fetching gene data...")

            for i, task in enumerate(tqdm(asyncio.as_completed(tasks), total=len(tasks), desc="Downloading Results", unit="gene")):
                try:
                    success = await task
                    if not success:
                        failed_genes.append(genes_to_fetch[i])
                except Exception as e:
                    print(f"{Fore.RED}Error processing gene {genes_to_fetch[i]}: {e}{Style.RESET_ALL}")
                    failed_genes.append(genes_to_fetch[i])

            if not failed_genes:
                print("All genes fetched successfully.")
                break

            retry_count += 1

        if failed_genes:
            print(f"\n{Fore.RED}Failed to fetch data for {len(failed_genes)} genes after {max_retries} retries.{Style.RESET_ALL}")
            print(f"Failed Genes: {failed_genes}")

def parse_taxa_list(file_path):
    """Parse a CSV/TSV file containing a list of taxa."""
    try:
        df = pd.read_csv(file_path)
        taxa_list = set(df['taxa'].str.strip())
        print(f"Parsed {len(taxa_list)} taxa from {file_path}")
        return taxa_list
    except Exception as e:
        print(f"Error parsing taxa list file: {e}")
        raise

def filter_uniprot_data(uniprot_tsv_file, taxa_list):
    """Load and filter UniProt data by the 'Organism' column, selecting one gene per organism."""
    print("Filtering UniProt data based on taxa list...")
    uniprot_df = pd.read_csv(uniprot_tsv_file)

    if 'Organism' not in uniprot_df.columns or 'Gene Names (primary)' not in uniprot_df.columns:
        raise ValueError("Required columns missing from the UniProt data.")

    uniprot_df['Normalized Gene'] = uniprot_df['Gene Names (primary)'].apply(normalize_gene_name)
    uniprot_df = uniprot_df.dropna(subset=['Normalized Gene'])
    uniprot_df = uniprot_df.drop_duplicates(subset=['Organism', 'Normalized Gene'])

    filtered_df = uniprot_df[uniprot_df['Organism'].isin(taxa_list)]
    if filtered_df.empty:
        print("No matching taxa found.")
        return None

    filtered_df.to_csv(FILTERED_TAXA_FILE, sep='\t', index=False)
    print(f"Filtered data saved to {FILTERED_TAXA_FILE}")
    return filtered_df

async def main(taxa=None, filter_only=False):
    """Main function to orchestrate the fetching and filtering process."""
    start_time = time.time()
    categories = ["J", "L", "C", "G", "O", "D", "M"]
    output_file = 'uniprot_data.csv'

    print(f"{Fore.YELLOW}** Running COG Fetcher: Gene Data Fetcher Script **{Style.RESET_ALL}")

    if filter_only:
        if not taxa:
            print(f"{Fore.RED}Error: --filter flag requires --taxa to be specified.{Style.RESET_ALL}")
            return
        print(f"{Fore.BLUE}Running in filter-only mode...{Style.RESET_ALL}")
        taxa_list = parse_taxa_list(taxa)
        filtered_data = filter_uniprot_data(output_file, taxa_list)
        if filtered_data is not None:
            print(f"Filtered data contains {len(filtered_data)} entries.")
        return

    print(f"{Fore.BLUE}Fetching Housekeeping Genes...{Style.RESET_ALL}")
    df = await fetch_all_categories(categories)

    print(f'Combined Results Length: {len(df)}')
    df['Organism'] = pd.to_numeric(df['Organism'], errors='coerce')
    filtered_df = df[df['Organism'] >= 2066]

    print(f'Filtered Results Length: {len(filtered_df)}')

    print(f"{Fore.BLUE}Fetching Gene Names for COG IDs...{Style.RESET_ALL}")
    async with aiohttp.ClientSession() as session:
        all_gene_names = await asyncio.gather(*[fetch_cog_gene_names(session, cog_id) for cog_id in filtered_df['COG']])

    all_gene_names = [gene for sublist in all_gene_names for gene in sublist]

    print(f'Normalized Gene Names: {len(all_gene_names)}')
    print(f'{all_gene_names[:10]}...')

    print(f"{Fore.BLUE}Fetching UniProt data...{Style.RESET_ALL}")
    await fetch_all_genes(all_gene_names, output_file)

    print(f"UniProt data written to {output_file}")

    if taxa and not filter_only:
        print(f"{Fore.BLUE}Running filtering after fetching...{Style.RESET_ALL}")
        taxa_list = parse_taxa_list(taxa)
        filtered_data = filter_uniprot_data(output_file, taxa_list)
        if filtered_data is not None:
            print(f"Filtered data contains {len(filtered_data)} entries.")

    print(f"{Fore.GREEN}Script Execution Completed!{Style.RESET_ALL}")
    print(f"{Fore.YELLOW}Total Execution Time: {(time.time() - start_time) / 60:.2f} minutes{Style.RESET_ALL}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch and filter UniProt data based on COG categories and taxa.")
    parser.add_argument("--taxa", type=str, help="Path to the CSV file containing the list of taxa.")
    parser.add_argument("--filter", action='store_true', help="Filter the UniProt data based on the provided taxa list.")
    args = parser.parse_args()

    asyncio.run(main(taxa=args.taxa, filter_only=args.filter))