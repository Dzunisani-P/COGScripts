# COG Database Construction and Curation

This repository contains a set of Python scripts for constructing and curating a non-redundant COG-based protein database from a list of COG gene identifiers. These scripts generate refined databases for the metaproteomic analysis of prokaryotic (bacterial) microbial taxa, for the purpose of microbiome related research.

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Script Descriptions](#script-descriptions)
- [License](#license)

---

## Overview

The purpose of this repository is to provide a pipeline for creating and curating a COG-based protein database. This workflow broadly includes the following steps:

1. **Initial Gene Query**: Query the UniProt API for genes listed in a CSV or TSV file.
2. **Filter Results**: Apply filtering steps to generate subsets based on taxonomic relevance.
3. **Refine Data**: Refine the data to reduce redundancy.
4. **Database Construction**: The final database is constructed using the protein sequences for each protein in the refined datasets for integration into downstream metaproteomic analysis.

---

## Requirements

Before running the scripts, ensure the following Python packages are installed:

- `aiohttp==3.11.10`
- `colorama==0.4.6`
- `pandas==2.2.3`
- `tqdm==4.66.2`

These packages can be installed using the provided `requirements.txt` file.

---

## Installation

1. **Clone the Repository**:  
   Clone this repository to your local machine using the following command:

   ```bash
   git clone https://github.com/Dzunisani-P/COGScripts.git

2. **Install Dependencies**:
    Navigate to the cloned repository's directory and install the dependencies using pip:

    ```bash
    cd COGScripts
    pip install -r requirements.txt

---

## Usage

### Running the Pipeline
The primary script, **cogCuration.py**, is used to automate the entire database construction and curation pipeline. It takes a CSV or TSV file containing gene indentifiers as input, see **COG_proteins.csv** for reference.

To use the pipeline, follow these steps:

1. **Prepare Input File**:
    Create or obtain a CSV or TSV file containing the COG genes you wish to query. The file should list gene identifiers that can be used to query the UniProt database.
2. **Run the Main Script**:
    Execute the cogCuration.py script with the input file as an argument:

    ```bash
    python cogCuration.py COG_proteins.csv

The script will:
- Query the uniprot protein database using async_geneQuery.py
- Filter results with queryFiltering.py
- Refine the data using dataRefinement.py
- Perform final sequence database construction using genePetch.py
- The curated database will be saved in the databases directory within the repository.

### Output
- A .FASTA containing a non-redundant COG-based protein sequence database will stored in the databases directory.
- A .CSV containing the RAW query data from the Uniprot protein database in the databases directory.
- Intermediate .CSV results and filtered datasets will be stored in the subsets directory.

---

## Script Descriptions

### cogCuration.py
Purpose: Orchestrates the entire COG database construction and curation process.
Usage: Takes an input file (CSV/TSV format) with gene data and runs the necessary scripts to query, filter, refine, and curate the COG-based protein database.
### async_geneQuery.py
Purpose: Queries the UniProt API asynchronously to gather gene data based on identifiers in the input file.
Usage: This script takes an input file containing gene identifiers and fetches the corresponding data from UniProt.
### queryFiltering.py
Purpose: Filters the queried gene data based on specific criteria such as taxonomic relevance.
Usage: This script processes the output from async_geneQuery.py to ensure the data is relevant and ready for further refinement.
### dataRefinement.py
Purpose: Refines the filtered gene data to remove redundancy and prepare it for the final curation step.
Usage: This script processes the filtered data and applies additional refinement steps for quality control.
### genePetch.py
Purpose: Finalizes the curation process by applying additional steps such as organizing the data and preparing it for integration into metaproteomic workflows.
Usage: This script performs the last step of curation and prepares the non-redundant COG-based protein database.

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.





