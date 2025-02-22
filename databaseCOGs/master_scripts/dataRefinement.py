import os
import pandas as pd
import subprocess
import time
from colorama import Fore, Style, init
import argparse

init(autoreset=True)

def get_csv_files(directory, keywords=None):
    """
    Retrieves a list of CSV files from the specified directory. Optionally filters the files based on provided keywords.
    
    Args:
        directory (str): The path to the directory containing the files.
        keywords (list, optional): A list of keywords to filter files by their names.
    
    Returns:
        list: A list of file names (strings) matching the CSV format and optional keyword filter.
    """
    return [file for file in os.listdir(directory) if file.endswith('.csv') and (not keywords or any(keyword in file for keyword in keywords))]

def collapse_duplicates(df, key_columns):
    """
    Removes duplicate rows from a DataFrame based on specified key columns. Keeps the first occurrence.
    
    Args:
        df (pandas.DataFrame): The DataFrame to remove duplicates from.
        key_columns (list): A list of column names used to identify duplicates.
    
    Returns:
        pandas.DataFrame: The DataFrame with duplicates removed.
    
    Raises:
        ValueError: If none of the specified key columns are found in the DataFrame.
    """
    for key in key_columns:
        if key in df.columns:
            return df.drop_duplicates(subset=[key], keep='first')
    raise ValueError("None of the key columns found in the DataFrame.")

def save_dataframes(dataframes, suffix, output_dir):
    """
    Saves a dictionary of DataFrames to CSV files, appending a suffix to the file names.
    
    Args:
        dataframes (dict): A dictionary where keys are file names and values are DataFrames to be saved.
        suffix (str): A suffix to append to the output file names.
        output_dir (str): The directory where the output CSV files will be saved.
    """
    for name, df in dataframes.items():
        file_path = os.path.join(output_dir, f'{name}_{suffix}.csv')
        if os.path.exists(file_path):
            os.remove(file_path)
        df.to_csv(file_path, index=False)

def prepare_column_for_overlap(df, column):
    """
    Prepares a column for overlap calculation by extracting the first two words of each value (if applicable).
    
    Args:
        df (pandas.DataFrame): The DataFrame containing the column.
        column (str): The column name to be processed.
    
    Returns:
        pandas.Series: A Series containing the first two words of each value in the specified column.
    
    Raises:
        ValueError: If the specified column is not found in the DataFrame.
    """
    if column in df.columns:
        return df[column].dropna().apply(lambda x: ' '.join(x.split()[:2]))
    raise ValueError(f"Column '{column}' not found.")

def calculate_overlaps(dataframes, column, output_dir, overlap_suffix='overlap_results.csv'):
    """
    Calculates the overlap between values in the specified column across multiple DataFrames.
    
    Args:
        dataframes (dict): A dictionary where keys are DataFrame names and values are the DataFrames to compare.
        column (str): The column name to check for overlaps.
        output_dir (str): The directory where the overlap results will be saved.
        overlap_suffix (str, optional): The suffix for the output file name. Defaults to 'overlap_results.csv'.
    
    Returns:
        pandas.DataFrame: A DataFrame containing the overlap results between the DataFrames.
    """
    column_sets = {name: set(prepare_column_for_overlap(df, column)) for name, df in dataframes.items()}
    overlap_results = [
        {'DF1': df1, 'DF2': df2, 'Overlap Count': len(set1.intersection(set2)), 'Overlap Items': ', '.join(set1.intersection(set2))}
        for df1, set1 in column_sets.items() for df2, set2 in column_sets.items() if df1 < df2
    ]
    overlap_df = pd.DataFrame(overlap_results)
    overlap_path = os.path.join(output_dir, overlap_suffix)
    if os.path.exists(overlap_path):
        os.remove(overlap_path)
    overlap_df.to_csv(overlap_path, index=False)
    print(f"{Fore.GREEN}Overlap results saved as {overlap_path}{Style.RESET_ALL}")
    return overlap_df

def print_filtered_summary(directory, suffix):
    """
    Prints a summary of line counts for filtered files in the specified directory with the given suffix.
    
    Args:
        directory (str): The directory to search for filtered files.
        suffix (str): The suffix to filter the files by.
    """
    filtered_files = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(suffix)]
    if filtered_files:
        command = f"wc -l {' '.join(filtered_files)}"
        subprocess.run(command, shell=True, capture_output=True, text=True)
    else:
        print(f"No files with suffix '{suffix}' found for summary.")

def collapse_and_calculate_overlaps(input_dir, output_dir):
    """
    Removes duplicates in the input CSV files and calculates overlaps based on the 'Gene/Org' and 'Gene/Id' columns.
    
    Args:
        input_dir (str): The directory containing input CSV files.
        output_dir (str): The directory where output files will be saved.
    
    Returns:
        tuple: A tuple containing:
            - pandas.DataFrame: The overlap results DataFrame.
            - int: The total number of records after collapsing duplicates.
    """
    csv_files = get_csv_files(input_dir, keywords=['gene_org', 'gene_id'])
    dataframes = {os.path.splitext(file)[0]: pd.read_csv(os.path.join(input_dir, file)) for file in csv_files}
    collapsed_dataframes, total_collapsed_records = {}, 0
    for name, df in dataframes.items():
        try:
            collapsed_df = collapse_duplicates(df, key_columns=['Gene/Org', 'Gene/Id'])
            collapsed_dataframes[name] = collapsed_df
            total_collapsed_records += len(collapsed_df)
            print(f"Collapsed: {name} ({len(collapsed_df)} records)")
        except ValueError as e:
            print(f"{Fore.RED}Skipping {name}: {e}{Style.RESET_ALL}")
    save_dataframes(collapsed_dataframes, suffix='collapsed', output_dir=output_dir)
    overlap_results = calculate_overlaps(collapsed_dataframes, column='Organism', output_dir=output_dir)
    return overlap_results, total_collapsed_records

def priority_filter_dataframes(input_dir, output_dir, collapsed_files, priority_order):
    """
    Filters the collapsed dataframes based on a priority order and removes overlaps.
    
    Args:
        input_dir (str): The directory containing input CSV files.
        output_dir (str): The directory where output files will be saved.
        collapsed_files (list): A list of collapsed CSV files to process.
        priority_order (list): A list of priority names to filter the DataFrames.
    
    Returns:
        tuple: A tuple containing:
            - dict: A dictionary of filtered DataFrames.
            - pandas.DataFrame: The overlap results after filtering.
            - int: The total number of records after filtering.
    """
    dataframes = {os.path.splitext(file)[0]: pd.read_csv(os.path.join(input_dir, file)) for file in collapsed_files}
    sorted_dataframes = sorted(dataframes.items(), key=lambda item: next((i for i, s in enumerate(priority_order) if s in item[0]), len(priority_order)))
    filtered_dataframes, seen_items, total_filtered_records = {}, set(), 0
    for name, df in sorted_dataframes:
        try:
            df['First_Two_Organism'] = prepare_column_for_overlap(df, 'Organism')
        except ValueError as e:
            print(f"{Fore.RED}Skipping {name}: {e}{Style.RESET_ALL}")
            continue
        if seen_items:
            df = df[~df['First_Two_Organism'].isin(seen_items)]
        
        output_path = os.path.join(output_dir, f'{name}_filtered.csv')
        if os.path.exists(output_path):
            os.remove(output_path)

        filtered_dataframes[name] = df
        seen_items.update(df['First_Two_Organism'])
        total_filtered_records += len(df)
        df.to_csv(output_path, index=False)
        print(f"Filtered: {name} ({len(df)} records)")

    print_filtered_summary(output_dir, suffix='_filtered.csv')
    overlap_df = calculate_overlaps(filtered_dataframes, column='First_Two_Organism', output_dir=output_dir, overlap_suffix='filtered_organism_overlap_results.csv')
    return filtered_dataframes, overlap_df, total_filtered_records


def main(input_dir, output_dir):
    start_time = time.time()
    print(f"{Fore.YELLOW}** Running Data Refinement Script **{Style.RESET_ALL}")

    os.makedirs(output_dir, exist_ok=True)

    print("Removing Duplicates..")
    overlap_results, total_collapsed_records = collapse_and_calculate_overlaps(input_dir, output_dir)

    print("Removing Overlaps....")
    collapsed_files = get_csv_files(output_dir, keywords=['collapsed'])
    priority_order = ['species_gene_id', 'genus_gene_org', 'sp_gene_org', 'no_rank_gene_org', 'other_gene_org']
    filtered_dataframes, filtered_overlap_df, total_filtered_records = priority_filter_dataframes(input_dir, output_dir, collapsed_files, priority_order)

    elapsed_time = round((time.time() - start_time) / 60, 2)
    print(f"Total records across collapsed DataFrames: {total_collapsed_records}")
    print(f"Total records across filtered DataFrames: {total_filtered_records}")
    print(f"{Fore.GREEN}Script Execution Completed {Style.RESET_ALL}")
    print(f"{Fore.CYAN}Elapsed Time: {elapsed_time} mins{Style.RESET_ALL}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Refine data by removing duplicates and overlaps.")
    parser.add_argument('--input_dir', type=str, default=os.getcwd(), help="Path to the directory containing input files.")
    parser.add_argument('--output_dir', type=str, default=os.getcwd(), help="Path to the directory for saving output files.")
    args = parser.parse_args()

    main(args.input_dir, args.output_dir)
