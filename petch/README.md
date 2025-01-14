# PETCH: Proteome Fetcher Script

PETCH is an asynchronous Python script to fetch, filter, and validate proteomes from UniProt for a given taxa list. The script downloads reference and other proteomes, searches for matches with the taxa list, and downloads the corresponding proteome FASTA files. It also validates the data by comparing proteome counts with FASTA entries.

## Features

- Fetch proteomes from UniProt using asynchronous requests.
- Filter proteomes based on taxonomic lineage.
- Match taxa with reference and other proteomes.
- Download FASTA files for matched taxa.
- Count FASTA entries per organism and validate proteome counts.

## Requirements

- Python 3.7+
- Required Python packages:
  - `aiohttp`
  - `pandas`
  - `tqdm`
  - `colorama`

Install the dependencies using:

```bash
pip install aiohttp pandas tqdm colorama

```
---

## Usage

Run the script with the following command:
   ```bash
   python petch_asyc.py <input_file> --filterout <keywords> --output_dir <output_directory>
   ```

**Arguments:**

input_file: CSV or TSV file containing a list of taxa.
--filterout: Optional list of keywords to filter out.
--output_dir: Directory to save output files (default: outputs).

---

## Workflow

1. **Download Proteomes:** Fetch reference and other proteomes from UniProt.
2. **Filter Taxa:** Load taxa list and filter based on user input.
3. **Find Matches:** Search for taxa matches in the downloaded proteomes.
4. **Fetch FASTA Files:** Download FASTA sequences for matched taxa.
5. **Count Entries:** Count FASTA entries per organism.
6. **Validate Data:** Compare proteome counts with FASTA entry counts.


### Output

matches.csv: Matched taxa and proteomes.
no_matches.csv: Taxa with no proteome match.
Petch.fasta: Downloaded FASTA sequences.
fasta_counts.csv: FASTA entry counts per organism.
valid_counts.csv: Validation results comparing proteome counts and FASTA entries.


### Example

```bash
python petch_asyc.py taxa_list.csv --filterout "Bacteria" --output_dir "proteomes_output"
```
---

## License

This project is licensed under the MIT License - see the LICENSE file for details.





