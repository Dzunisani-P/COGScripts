import argparse
import aiohttp
import asyncio
import os
import time
from tqdm import tqdm
from colorama import Fore, Style, init
import pandas as pd

init(autoreset=True)

async def fetch_fasta(session, protein_entry):
    """
    Fetches the FASTA sequence for a given protein entry from UniProt.
    
    Args:
        session (aiohttp.ClientSession): The session used for making HTTP requests.
        protein_entry (str): The UniProt protein entry identifier.

    Returns:
        str or None: The FASTA sequence if the request is successful, otherwise None.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{protein_entry}.fasta"
    try:
        async with session.get(url) as response:
            return await response.text() if response.status == 200 else None
    except Exception as e:
        print(f"{Fore.RED}Error fetching {protein_entry}: {e}{Style.RESET_ALL}")
        return None

async def process_batch(batch, session, progress_bar):
    """
    Processes a batch of protein entries by fetching their corresponding FASTA sequences.
    
    Args:
        batch (list): A list of protein entry identifiers to be processed.
        session (aiohttp.ClientSession): The session used for making HTTP requests.
        progress_bar (tqdm): The progress bar for tracking the download progress.

    Returns:
        filter: A filter object containing valid FASTA sequences (non-None results).
    """
    tasks = [fetch_fasta(session, entry) for entry in batch]
    results = await asyncio.gather(*tasks)
    progress_bar.update(len(batch))
    return filter(None, results)

async def download_fasta(entries, output_file, batch_size=50000):
    """
    Downloads FASTA sequences for a list of protein entries in batches and saves them to a file.
    
    Args:
        entries (list): A list of protein entry identifiers.
        output_file (str): The file path where the FASTA sequences will be saved.
        batch_size (int, optional): The number of entries to download in each batch.
    Returns:
        None
    """
    if os.path.exists(output_file):
        os.remove(output_file)  # Clear the file before starting
    async with aiohttp.ClientSession() as session:
        with tqdm(total=len(entries), desc="Downloading sequences", unit="entry") as progress_bar:
            for i in range(0, len(entries), batch_size):
                batch = entries[i:i+batch_size]
                fasta_sequences = await process_batch(batch, session, progress_bar)
                with open(output_file, 'a') as f:
                    f.write("\n".join(fasta_sequences) + "\n")

def get_entries(directory):
    """
    Retrieves the unique protein entry identifiers from CSV files in the specified directory.
    
    Args:
        directory (str): The directory containing the filtered CSV files.
    
    Returns:
        list: A list of unique protein entry identifiers found in the CSV files.
    """
    entries = []
    for file in os.listdir(directory):
        if file.endswith('filtered.csv'):
            try:
                df = pd.read_csv(os.path.join(directory, file))
                if 'Entry' in df.columns:
                    entries.extend(df['Entry'].dropna())
            except Exception as e:
                print(f"{Fore.RED}Error processing {file}: {e}{Style.RESET_ALL}")
    return list(set(entries))

def validate_fasta(output_file, expected_count):
    """
    Validates that the downloaded FASTA file contains the expected number of entries.
    
    Args:
        output_file (str): The path to the FASTA file to be validated.
        expected_count (int): The expected number of protein entries in the FASTA file.
    
    Returns:
        None
    """
    if not os.path.exists(output_file):
        print(f"{Fore.RED}FASTA file not found!{Style.RESET_ALL}")
        return
    with open(output_file, 'r') as f:
        fasta_entries = sum(1 for line in f if line.startswith(">"))
    if fasta_entries == expected_count:
        print(f"{Fore.GREEN}Validation Passed: {fasta_entries} entries in the FASTA.{Style.RESET_ALL}")
    else:
        print(f"{Fore.RED}Validation Failed: {fasta_entries} entries found, expected {expected_count}.{Style.RESET_ALL}")

def main(input_directory, output_directory):
    start_time = time.time()
    print(f"{Fore.YELLOW}** Running gene FASTA Download Script **{Style.RESET_ALL}")

    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)
    output_fasta = os.path.join(output_directory, "COG_sequence.fasta")

    entries = get_entries(input_directory)
    print(f"Total entries: {len(entries)}")

    print(f"Downloading FASTA sequences to: {output_fasta}")
    asyncio.run(download_fasta(entries, output_fasta, batch_size=50000))

    validate_fasta(output_fasta, len(entries))
    elapsed_time = round((time.time() - start_time) / 60, 2)

    print(f"{Fore.GREEN}Script Execution Completed{Style.RESET_ALL}")
    print(f"{Fore.CYAN}Elapsed Time: {elapsed_time} mins{Style.RESET_ALL}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download FASTA sequences for protein entries.")
    parser.add_argument("--input_dir", help="Directory containing the filtred CSV files.")
    parser.add_argument("--output_dir", help="Directory to save the downloaded FASTA sequences.")
    args = parser.parse_args()

    main(args.input_dir, args.output_dir)
