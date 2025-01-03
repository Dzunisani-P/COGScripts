import aiohttp
import asyncio
import gzip
from io import StringIO
import pandas as pd
from aiohttp import ClientSession
from tqdm.asyncio import tqdm
import argparse
import os
import time
from colorama import Fore, Style

def read_cog_genes(input_file):
    """
    Reads a CSV/TSV file containing COG gene identifiers and returns a list of genes.

    Args:
        input_file (str): Path to the input CSV/TSV file.

    Returns:
        list: A list of gene names from the file.

    Raises:
        ValueError: If the file does not contain a 'genes' column or if the file format is unsupported.
    """

    if input_file.endswith(('.csv', '.tsv')):
        sep = ',' if input_file.endswith('.csv') else '\t'
        cogs_df = pd.read_csv(input_file, sep=sep)

        gene_column = cogs_df.columns[cogs_df.columns.str.lower() == 'genes']
        if gene_column.empty:
            raise ValueError("File does not contain a 'genes' column.")
        
        genes_column_name = gene_column[0]
        print(f"COG Genes Found: {len(cogs_df)}")
        return cogs_df[genes_column_name].tolist()
    
    raise ValueError("Unsupported file format. Please provide a CSV or TSV file.")


async def fetch_gene_data(session, gene, output_file, retries=3, batch_size=100):
    """
    Fetches gene data from UniProt's API for a single gene with pagination.

    Args:
        session (ClientSession): The aiohttp session used to make requests.
        gene (str): The gene to query data for.
        output_file (str): Path to the output CSV file to save the data.
        retries (int, optional): Number of retry attempts in case of failure.
        batch_size (int, optional): Number of records to write to the output file in each batch.

    Returns:
        bool: True if data fetching was successful, False otherwise.
    """

    query_url = f"https://rest.uniprot.org/uniprotkb/search?compressed=true&fields=accession%2Cxref_proteomes%2Cid%2Cgene_names%2Corganism_name%2Corganism_id%2Clineage_ids&format=tsv&query=gene:{gene}&size=500"
    
    batch_data = []
    try:
        while query_url:
            async with session.get(query_url) as response:
                response.raise_for_status()
                decompressed_data = gzip.decompress(await response.read()).decode('utf-8')
                df = pd.read_csv(StringIO(decompressed_data), sep='\t')
                batch_data.append(df)

                if len(batch_data) >= batch_size:
                    df_batch = pd.concat(batch_data, ignore_index=True)
                    df_batch.to_csv(output_file, mode='a', header=not os.path.exists(output_file), index=False)
                    batch_data = [] 
                del df

                next_page_url = response.headers.get("link")
                query_url = next_page_url.split(";")[0].strip("<>") if next_page_url and 'rel="next"' in next_page_url else None
        
        if batch_data:
            df_batch = pd.concat(batch_data, ignore_index=True)
            df_batch.to_csv(output_file, mode='a', header=not os.path.exists(output_file), index=False)
        
        return True

    except Exception as e:
        if retries > 0:
            await asyncio.sleep(2 ** (3 - retries))  
            return await fetch_gene_data(session, gene, output_file, retries - 1)
        print(f"Failed to fetch data for gene {gene}: {e}")
        return False


async def fetch_all_genes(cog_genes, output_file):
    """
    Fetches data for all COG genes concurrently using asynchronous requests.

    Args:
        cog_genes (list): A list of COG gene names.
        output_file (str): Path to the output CSV file to save the data.

    Returns:
        None: Prints status updates and information about failed genes to the console.
    """

    if os.path.exists(output_file):
        os.remove(output_file)
        print(f"Existing output file {output_file} removed.")

    async with ClientSession() as session:
        tasks = [fetch_gene_data(session, gene, output_file) for gene in cog_genes]
        print("Fetching gene data...")

        failed_genes = []
        for i, task in enumerate(tqdm(asyncio.as_completed(tasks), total=len(tasks), desc="Downloading Results", unit="gene")):
            try:
                success = await task
                if not success:
                    failed_genes.append(cog_genes[i])
            except Exception as e:
                print(f"{Fore.RED}Error processing gene {cog_genes[i]}: {e}{Style.RESET_ALL}")
                failed_genes.append(cog_genes[i])

        if failed_genes:
            print(f"\n{Fore.RED}Failed to fetch data for {len(failed_genes)} genes.{Style.RESET_ALL}")
            print(f"Failed Genes: {failed_genes}")



# Main function
async def main(input_file, output_dir):
    start_time = time.time()

    print(f"{Fore.YELLOW}** Running COG Gene Querying Script **{Style.RESET_ALL}")
    cog_genes = read_cog_genes(input_file)

    if not output_dir:
        output_dir = os.getcwd()
        
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, os.path.basename(input_file).replace(os.path.splitext(input_file)[1], '_gene_queries.csv'))

    await fetch_all_genes(cog_genes, output_file)

    print(f"{Fore.GREEN}Script Execution Completed.{Style.RESET_ALL}")
    print(f"{Fore.CYAN}Elapsed Time: {(time.time() - start_time) / 60:.2f} mins{Style.RESET_ALL}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch UniProt query data for genes from a CSV/TSV file.")
    parser.add_argument("input_file", help="Input CSV or TSV file containing genes.")
    parser.add_argument("--output_dir", help="Directory to save the output CSV file. Defaults to the current working directory.", default=None)
    args = parser.parse_args()

    asyncio.run(main(args.input_file, args.output_dir))

