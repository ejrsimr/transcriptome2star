#!/usr/bin/env python:w!
# transcriptome2star.py
# author: Eric Ross (ejr@stowers.org)
# Developed in order to use transcriptomes as single cell references in collaboration with Cynthia Chen
# convert the FASTA file of a transcriptome into a FASTA and GTF file suitable for creating a STAR reference
# usage: python transcriptome2star.py <fasta_in> <gtf_out> <fasta_out>

# import modules
import sys

def main():
    
    # get arguments
    if len(sys.argv) != 4:
        print("usage: python transcriptome2star.py <fasta_in> <gtf_out> <fasta_out>")
        return  # Exit immediately if the number of arguments is incorrect

    fasta_in = sys.argv[1]
    gtf_out = sys.argv[2]
    fasta_out = sys.argv[3]
    
    # check if files exist
    # check if input file exists
    try:
        with open(fasta_in, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: {fasta_in} not found")
        exit(1)
    # check if output files already exist
    try:
        with open(gtf_out, 'r') as f:
            print(f"Error: {gtf_out} already exists")
            exit(1)
    except FileNotFoundError:
        pass
    try:
        with open(fasta_out, 'r') as f:
            print(f"Error: {fasta_out} already exists")
            exit(1)
    except FileNotFoundError:
        pass
    
    # read in FASTA file
    fasta = read_fasta(fasta_in)
    
    # OUTPUT GTF and FASTA files
    print_gtf(fasta, gtf_out)
    print_fasta(fasta, fasta_out)
    
    exit(0)

def read_fasta(filename):
    """
    read in the fasta file
    """
    fasta = {}
    header = None
    # load fasta into dictionary
    with open(filename, 'r') as input_file:
        for line in input_file:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                header = header.split(' ')[0]
            else:
                fasta.setdefault(header, "")
                fasta[header] += line
    return fasta

def print_gtf(fasta, gtf_out):
    """
    print the gtf file
    """
    # print fasta as gtf: ordered by sequence length
    with open(gtf_out, 'w') as output_file:
        for header, sequence in sorted(fasta.items(), key=lambda x: len(x[1]), reverse=True): 
            seq_length = len(sequence)
            output_file.write(f"{header}.ref\tejr\ttranscript\t1\t{seq_length}\t.\t+\t.\tgene_id \"{header}\"; transcript_id \"{header}\"; gene_name \"{header}\";\n")
            output_file.write(f"{header}.ref\tejr\texon\t1\t{seq_length}\t.\t+\t.\tgene_id \"{header}\"; transcript_id \"{header}\"; gene_name \"{header}\"; exon_number \"1\"; exon_id \"{header}.1\";\n")

def print_fasta(fasta, fasta_out):
    """
    print the fasta file
    """
    # print the fasta file
    with open(fasta_out, 'w') as output_file:
        for header in sorted(fasta, key=lambda k: len(fasta[k]), reverse=True):
            # print the sequence
            sequence = fasta[header]
            output_file.write(f">{header}.ref\n")
            sequence = [sequence[i:i+60] for i in range(0, len(sequence), 60)]
            output_file.write('\n'.join(sequence) + '\n')      

# run main
if __name__ == "__main__":
    main()