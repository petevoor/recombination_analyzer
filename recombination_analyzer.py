#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 11:57:11 2024

@author: petevoor
"""

import argparse
from Bio import SeqIO
import pandas as pd
import os
from tqdm import tqdm
import subprocess
import tempfile

attL = "TCGGCCGGCTTGTCGACGACGGCGTCCTCAGTGGTGTACGGTACAAACC"
sub_seq_str = "gtcgagggtgaagtacttgctgacttccttgaggaacacatgatgcgtcctacggttgctgctacgcatatcattgagatgtctgtgggaggagttgatgtgtactctgaggacgatgagggttacggtacgtctttcattgagtggtgatttatgcattaggactgcatagggatgcactatagaccacggatggtcagttctttaagttactgaaaagacacgataaattaatacgactcactatagggagaggagggacgaaaggttactatatagatactgaatgaatacttatagagtgcataaagtatgcataatggtgtacctagagtgacctctaagaatggtgattatattgtattagtatcaccttaacttaagttcctatagggtcctttaaaatataccataaaaatctgagtgactatctcacagtgtacggacctaaagttcccccatagggggtacctaaagcccagccaatcacctaaagtcaaccttcggttgaccttgagggttccctaagggttggggatgacccttgggtttgtctttgggtgttaccttgagtgtctctctgtgtccctatctgttacagtctcctaaagtatcctcctaaagtcacctcctaacgctcgaggccctat"

def create_blast_database(reference_path):
    cmd = ["makeblastdb", "-in", reference_path, "-dbtype", "nucl"]
    subprocess.run(cmd, check=True)

def align_with_blast(read_seq, reference_path):
    cmd = ["blastn", "-db", reference_path, "-query", "-", "-outfmt", "6 qseqid pident length"]
    input_seq = f">query\n{str(read_seq)}"
    result = subprocess.run(cmd, input=input_seq, text=True, capture_output=True, check=False)
    lines = result.stdout.strip().split("\n")
    if lines and len(lines[0].split()) == 3:  # Ensure there's an alignment and it contains 3 fields
        _, pident, length = lines[0].split()
        return float(pident), int(length)
    else:
        return 0, 0
    
def align_subsequence(read_seq, sub_seq_str):
    # Create a temporary file to hold the subsequence
    with tempfile.NamedTemporaryFile(mode='w+', delete=True) as temp_file:
        temp_file.write(f">subsequence\n{sub_seq_str}\n")
        temp_file.flush()  # Ensure data is written to disk

        cmd = ["blastn", "-query", "-", "-subject", temp_file.name, "-outfmt", "6 qseqid pident length"]
        input_seq = f">query\n{str(read_seq)}"
        result = subprocess.run(cmd, input=input_seq, text=True, capture_output=True, check=False)
        
        lines = result.stdout.strip().split("\n")
        if lines and len(lines[0].split()) == 3:  # Ensure there's an alignment and it contains 3 fields
            _, pident, _ = lines[0].split()
            return float(pident)
        else:
            return 0.0

def process_fastq(fastq_file, reference_path, sub_seq_str):
    total_reads = 0
    matching_reads = 0
    sub_seq_reads = 0
    diagnostic_data = []

    for record in tqdm(SeqIO.parse(fastq_file, "fastq"), desc="Processing reads", unit="read", leave=False, position=1):
        total_reads += 1
        pident, aligned_length = align_with_blast(record.seq, reference_path)

        if aligned_length >= 0.7 * len(record.seq) and pident > 80:
            sub_pident = align_subsequence(record.seq, sub_seq_str)
            
            if sub_pident >= 75:
                matching_reads += 1
                sub_seq_reads += 1
                
            else:
                attL_pident = align_subsequence(record.seq, attL)
                
                if attL_pident >= 99:
                    matching_reads += 1
                    

            diagnostic_data.append([record.id, len(record.seq), pident, sub_pident])

    return total_reads, matching_reads, sub_seq_reads, diagnostic_data


def main(reference_file, fastq_dir, output_file_name):
    
    create_blast_database(reference_file)

    fastq_files = [os.path.join(fastq_dir, f) for f in os.listdir(fastq_dir) if f.endswith('.fastq')]

    results = []
    all_diagnostic_data = []

    for fastq_file in tqdm(fastq_files, desc="Processing files", unit="file", leave=False, position=0):
        total, matches, subs, diagnostic_data = process_fastq(fastq_file, reference_file, sub_seq_str)
        percent_sub = (subs / matches) * 100 if matches != 0 else 0
        results.append([fastq_file, total, matches, subs, percent_sub])
        all_diagnostic_data.extend(diagnostic_data)

    df = pd.DataFrame(results, columns=["File Name", "Total Reads", "Matching Reads", "Sub-sequence Reads", "Percent Sub-sequence"])
    
    output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), output_file_name)
    df.to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process FASTQ files against a reference sequence.')
    parser.add_argument('reference_path', type=str, help='Path to the reference sequence')
    parser.add_argument('fastq_dir', type=str, help='Directory containing FASTQ files')
    parser.add_argument('output_file_name', type=str, help='Name of the output CSV file')

    args = parser.parse_args()

    main(args.reference_path, args.fastq_dir, args.output_file_name)
