
# Recombination Analyzer

This Python script processes FASTQ reads against a reference nucleotide sequence to identify and analyze whether it has been recombined or not. It leverages BLAST for sequence alignment and provides a comprehensive output in CSV format.

## Overview

- Creates a BLAST database from a provided reference sequence.
- Aligns each read from FASTQ files against the reference sequence.
- Identifies reads with at least 75% identity to the reference.
- If those reads contain the attP/attB flanked subsequence, they are counted as not recombined
- If those reads do not contain the subsequence and do contain an attL site, they are counted as recombined 
- Outputs a detailed summary including total reads, matching reads, subsequence-specific reads, and percentages.

## Dependencies

- Python 3.x
- Biopython
- Pandas
- tqdm
- BLAST

## Installation

1. Ensure Python 3.x is installed.
2. Install the required Python libraries:
   ```bash
   pip install biopython pandas tqdm
   ```
3. Install BLAST command-line tools from [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

## Usage

1. Prepare a reference nucleotide sequence file in .fasta format.
2. Place FASTQ files in a designated directory.
3. Run the script using the command line:
   ```bash
   python fastq_sequence_analyzer.py <reference_path> <fastq_dir> <output_file_name>
   ```
   Where:
   - `<reference_path>` is the path to the reference sequence file.
   - `<fastq_dir>` is the directory containing FASTQ files.
   - `<output_file_name>` is the desired name for the output CSV file.

## Output

The script generates a CSV file containing:
- File Name: Name of the processed FASTQ file.
- Total Reads: Total number of reads in the FASTQ file.
- Matching Reads: Number of reads matching the reference sequence.
- Sub-sequence Reads: Number of reads containing a specific subsequence.
- Percent Sub-sequence: Percentage of matching reads that contain the subsequence.

## Notes

- The script assumes the presence of a predefined subsequence `sub_seq_str` and `attL` sequence for alignment checks.
- Ensure that the BLAST tools are correctly installed and accessible from the command line.
