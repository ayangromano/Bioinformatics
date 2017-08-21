#!/usr/local/bin/python

from Bio import SeqIO
import argparse

def sort_scaffold(scaffold_file, output_name):
    """
    This function will sort scaffold fasta files by length from longest to shortest. 
    """
    records = list(SeqIO.parse(scaffold_file, 'fasta'))
    records.sort(cmp=lambda x, y: cmp(len(y), len(x)))
    SeqIO.write(records, output_name, 'fasta')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Enter input file name', required=True)
    parser.add_argument('-o', '--output', help='Enter output file name', required=True)
    args = parser.parse_args()

    sort_scaffold(args.input, args.output)

if __name__=="__main__":
    main()
