#!/usr/bin/env python3

import csv
import os

def schema_to_fasta(input_tsv, schema_name):
    input_tsv_dialect = csv.Sniffer().sniff(input_tsv.readline())
    input_tsv.seek(0)
    input_tsv_reader = csv.DictReader(input_tsv, dialect=input_tsv_dialect)
    output_fasta_path = f"{os.path.splitext(input_tsv.name)[0]}.fasta"
    with open(output_fasta_path, "w") as output_fasta:
        for row in input_tsv_reader:
            print(f">{schema_name}_{row['Amplicon']}_{row['Primer']}", file=output_fasta)
            print(row["Sequence"], file=output_fasta)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Convert a TSV file of primers into a FASTA file and two BED files (one for primers, one for amplicons).")
    parser.add_argument("input_tsv", type=argparse.FileType("r"), help="TSV file of primers.")
    parser.add_argument("--schema_name", type=str, help="Schema name (for FASTA output)")
    # parser.add_argument("--reference", type=str, help="Reference name (for BED outputs)")

    args = parser.parse_args()
    schema_to_fasta(args.input_tsv, args.schema_name)
