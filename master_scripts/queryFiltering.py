import pandas as pd
import re
import os
import argparse
import time
from colorama import Fore, Style

# Suppress the SettingWithCopyWarning globally
pd.set_option('mode.chained_assignment', None)

def filter_queries(input_csv):
    """
    Filters the input uniprot.csv to retain only bacterial taxa and cleans the 'Organism' column.

    Args:
        input_csv (str): Path to the input CSV file containing taxonomic data.

    Returns:
        pd.DataFrame: A DataFrame with only bacterial taxa and the 'Organism' column cleaned.
    """
    df = pd.read_csv(input_csv)
    print(f"Original data: {len(df)} records")
    print("Removing non-proteome and non-bacterial records..")
    df = df[df['Proteomes'].notna()]
    df = df[df['Taxonomic lineage (Ids)'].str.contains(r'2 \(superkingdom\)', na=False)]
    df['Organism'] = df['Organism'].apply(lambda x: re.sub(r'\(.*?\)', '', x))
    return df


def separate_organisms_with_sp(df):
    """
    Separates organisms whose names end with 'sp.' or 'sp' into a separate DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame containing the organism data.

    Returns:
        tuple: A tuple containing two DataFrames:
            - sp_rows: DataFrame with organisms ending in 'sp.' or 'sp'.
            - remaining_rows: DataFrame with all other organisms.
    """
    sp_rows = df[df['Organism'].str.strip().str.endswith(('sp.', 'sp'), na=False)]
    remaining_rows = df[~df['Organism'].str.strip().str.endswith(('sp.', 'sp'), na=False)]
    return sp_rows, remaining_rows


def separate_by_highest_taxonomic(df):
    """
    Classifies rows based on the last entry in the 'Taxonomic lineage (Ids)' column.

    Args:
        df (pd.DataFrame): The DataFrame containing taxonomic data.

    Returns:
        tuple: A tuple containing four DataFrames:
            - genus: Rows classified as genus.
            - species: Rows classified as species.
            - no_rank: Rows classified as 'no rank'.
            - other: Rows classified as other taxonomic levels.
    """
    genus, species, no_rank, other = [], [], [], []

    for _, row in df.iterrows():
        last_entry = row['Taxonomic lineage (Ids)'].split(',')[-1].strip().lower()
        if 'genus' in last_entry:
            genus.append(row)
        elif 'species' in last_entry:
            species.append(row)
        elif 'no rank' in last_entry:
            no_rank.append(row)
        else:
            other.append(row)

    return pd.DataFrame(genus), pd.DataFrame(species), pd.DataFrame(no_rank), pd.DataFrame(other)

# OPTIONAL validation function for checking dfs.
def get_taxonomic_level_counts(df):
    """
    Counts unique taxonomic levels from the last entry in the 'Taxonomic lineage (Ids)' column.

    Args:
        df (pd.DataFrame): The DataFrame containing taxonomic data.

    Returns:
        tuple: A tuple containing:
            - level_counts: A Series of counts of unique taxonomic levels.
            - total_count: The total number of unique taxonomic levels.
    """
    level_counts = df['Taxonomic lineage (Ids)'].str.split(',').str[-1].str.strip().str.extract(r'\((.*?)\)')[0].value_counts()
    total_count = level_counts.sum()
    return level_counts, total_count


def add_gene_org_column(df, df_name, output_dir):
    """
    Adds a 'Gene/Org' column combining the first gene name and the first two organism names.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        df_name (str): The name of the DataFrame to use for the output file.
        output_dir (str): The directory where the output CSV file will be saved.

    Returns:
        None: The function saves the DataFrame to a CSV file with the name <df_name>_gene_org.csv.
    """

    df.loc[:, 'Gene/Org'] = df.apply(
        lambda row: f"{row['Gene Names'].split()[0] if pd.notna(row['Gene Names']) else ''} "
                    f"{' '.join(row['Organism'].split()[:2]) if pd.notna(row['Organism']) else ''}",
        axis=1
    )
    output_path = os.path.join(output_dir, f"{df_name}_gene_org.csv")
    if os.path.exists(output_path):
        os.remove(output_path)
    df.to_csv(output_path, index=False)


def add_gene_id_column(df, df_name, output_dir):
    """
    Adds a 'Gene/Id' column combining the first gene name and the number before '(species)' in the 'Taxonomic lineage (Ids)' column.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        df_name (str): The name of the DataFrame to use for the output file.
        output_dir (str): The directory where the output CSV file will be saved.

    Returns:
        None: The function saves the DataFrame to a CSV file with the name <df_name>_gene_org.csv.
    """
    df.loc[:, 'Gene/Id'] = df.apply(
        lambda row: f"{row['Gene Names'].split()[0] if pd.notna(row['Gene Names']) else ''} "
                    f"{next((re.search(r'(\d+)', part).group(1) for part in row['Taxonomic lineage (Ids)'].split(',')
                           if '(species)' in part.lower()), '')}",
        axis=1
    )
    output_path = os.path.join(output_dir, f"{df_name}_gene_id.csv")
    if os.path.exists(output_path):
        os.remove(output_path)
    df.to_csv(output_path, index=False)



def main(input_file, output_dir):
    start_time = time.time()

    os.makedirs(output_dir, exist_ok=True)

    print(f"{Fore.YELLOW}** Running Query Processing Script **{Style.RESET_ALL}")
    # Filter and preprocess the data
    filtered_df = filter_queries(input_file)
    print(f"Filtered data: {len(filtered_df)} records")

    # Separate the organisms ending with 'sp.' or 'sp' from the rest
    sp_rows, remaining_rows = separate_organisms_with_sp(filtered_df)
    print(f"Sp. data: {len(sp_rows)} records")
    
    # Separate the remaining data into different taxonomic levels
    print("Separating remaining data..")
    genus_df, species_df, no_rank_df, other_df = separate_by_highest_taxonomic(remaining_rows)
    print(f"Genus ({len(genus_df)}), Species ({len(species_df)}), No Rank ({len(no_rank_df)}), Other ({len(other_df)})")
    print(f'Combined data: {len(genus_df) + len(species_df) + len(no_rank_df) + len(other_df)}')

    # Apply the relevant functions to each DataFrame
    print(f"Adding Gene/Org/ID columns..")
    genus_df = add_gene_org_column(genus_df, 'genus', output_dir)
    sp_rows = add_gene_org_column(sp_rows, 'sp', output_dir)

    # Process species data separately using add_gene_id_column
    species_df = add_gene_id_column(species_df, 'species', output_dir)
    other_df = add_gene_id_column(other_df, 'other', output_dir)
    no_rank_df = add_gene_id_column(no_rank_df, 'no_rank', output_dir)

    print(f"{Fore.GREEN}Script Execution Completed.{Style.RESET_ALL}")
    print(f"{Fore.CYAN}Elapsed Time: {(time.time() - start_time) / 60:.2f} mins{Style.RESET_ALL}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process bacterial taxonomic data.")
    parser.add_argument("input_file", help="Input CSV file containing bacterial data.")
    parser.add_argument("--output_dir", default=os.getcwd(), help="Directory to save output files.")
    args = parser.parse_args()

    main(args.input_file, args.output_dir)
