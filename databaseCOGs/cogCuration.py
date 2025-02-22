import os
import shutil
import time
import glob
import subprocess
from colorama import Fore, Style
import argparse

def run_script(script_name, *args):
    """
    Helper function to run a Python script as an external process.
    
    Args:
        script_name (str): The name of the script to run.
        *args: Additional arguments to pass to the script.
    
    Returns:
        bool: True if the script executed successfully, False otherwise.
    """
    try:
        subprocess.run(["python3", os.path.join(os.getcwd(), "master_scripts", script_name)] + list(args), check=True)
    except subprocess.CalledProcessError as e:
        print(f"{Fore.RED}Error in {script_name}: {e}{Style.RESET_ALL}")
        return False
    return True

def main(input_file):
    """
    Main function to construct and curate a non-redundant COG-based protein database by running a series of scripts.
    
    Args:
        input_file (str): The path to the input CSV or TSV file containing gene data.
    
    Returns:
        None
    """
    start_time = time.time()
    cwd = os.getcwd()

    print(f"{Fore.MAGENTA}COG Database Construction and Curation:\n{Style.RESET_ALL}")

    db_dir = os.path.join(cwd, 'databases')
    if os.path.exists(db_dir):
            shutil.rmtree(db_dir)
    os.makedirs(db_dir, exist_ok=True)
    input_file_path = os.path.join(cwd, input_file)

    # Run async_geneQuery.py script
    if not run_script("async_geneQuery.py", input_file_path, "--output_dir", db_dir):
        return

    # Check if the output file exists
    output_file = os.path.join(db_dir, "*_gene_queries.csv")
    matching_file = glob.glob(output_file)
    if not matching_file:
        print(f"{Fore.RED}No matching files found in {db_dir}.{Style.RESET_ALL}")
        return

    subsets_dir = os.path.join(db_dir, 'subsets')
    os.makedirs(subsets_dir, exist_ok=True)

    # Run the other scripts
    if not run_script("queryFiltering.py", matching_file[0], "--output_dir", subsets_dir):
        return

    if not run_script("dataRefinement.py", "--input_dir", subsets_dir, "--output_dir", subsets_dir):
        return

    if not run_script("genePetch.py", "--input_dir", subsets_dir, "--output_dir", db_dir):
        return

    print(f"{Fore.MAGENTA}Master Script Completed in {(time.time() - start_time) / 60:.2f} minutes.{Style.RESET_ALL}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct and curate a non-redundant COG-based protein database.")
    parser.add_argument("input_file", help="Input CSV or TSV file containing genes.")
    args = parser.parse_args()
    
    main(args.input_file)
