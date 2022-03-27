#!/usr/bin/env python3

import csv
import os
from collections import namedtuple
from dataclasses import dataclass
from enum import Enum


def schema_to_fasta(input_tsv, schema_name):
    input_tsv_dialect = csv.Sniffer().sniff(input_tsv.readline())
    input_tsv.seek(0)
    input_tsv_reader = csv.DictReader(input_tsv, dialect=input_tsv_dialect)
    output_fasta_path = f"{os.path.splitext(input_tsv.name)[0]}.fasta"
    with open(output_fasta_path, "w") as output_fasta:
        for row in input_tsv_reader:
            print(f">{schema_name}_{row['Amplicon']}_{row['Primer']}", file=output_fasta)
            print(row["Sequence"], file=output_fasta)


class Sense(Enum):
    POSITIVE = "+"
    NEGATIVE = "-"


@dataclass
class Primer:
    direction: Sense
    start: int
    end: int
    pool: int
    name_alias: str


Amplicon = namedtuple("Amplicon", "start end pool")


def schema_to_beds(input_tsv, schema_name, reference_name):
    input_tsv.seek(0)
    input_tsv_dialect = csv.Sniffer().sniff(input_tsv.readline())
    input_tsv.seek(0)
    input_tsv_reader = csv.DictReader(input_tsv, dialect=input_tsv_dialect)
    output_primer_bed_path = f"{os.path.splitext(input_tsv.name)[0]}.primer.bed"
    output_amplicon_bed_path = f"{os.path.splitext(input_tsv.name)[0]}.amplicon.bed"

    amplicon_primers = dict()
    for row in input_tsv_reader:
        # breakpoint()
        match row["Sense"]:
            case "+":
                this_sense = Sense.POSITIVE
            case "-":
                this_sense = Sense.NEGATIVE
            case _:
                raise ValueError(f"invalid sense for amplicon {row['Amplcion']}: {row['Sense']}")
        if row["Amplicon"] not in amplicon_primers:
            amplicon_primers[row["Amplicon"]] = list()
        amplicon_primers[row["Amplicon"]].append(Primer(this_sense, int(row["Start"]), int(row["End"]), int(row["Pool"]), row["Primer"]))

    with open(output_primer_bed_path, "w") as output_primer_bed:
        for amplicon, primers in amplicon_primers.items():
            for primer in primers:
                primer_name = f"{schema_name}_{amplicon}_{primer.name_alias}"
                primer_bed_row = [
                    reference_name,
                    str(primer.start),
                    str(primer.end),
                    primer_name,
                    str(primer.pool),
                    primer.direction.value
                ]
                print("\t".join(primer_bed_row), file=output_primer_bed)

    amplicons = dict()
    for amplicon, primers in amplicon_primers.items():
        positive_sense_primers = [primer for primer in primers if primer.direction == Sense.POSITIVE]
        positive_sense_starts = [primer.start for primer in positive_sense_primers]
        max_start = max(positive_sense_starts)
        positive_sense_ends = [primer.end for primer in positive_sense_primers]
        max_end = max(positive_sense_ends)
        amp_start = max(max_start, max_end)

        negative_sense_primers = [primer for primer in primers if primer.direction == Sense.NEGATIVE]
        negative_sense_starts = [primer.start for primer in negative_sense_primers]
        min_start = min(negative_sense_starts)
        negative_sense_ends = [primer.end for primer in negative_sense_primers]
        min_end = min(negative_sense_ends)
        amp_end = min(min_start, min_end)

        amplicons[amplicon] = Amplicon(amp_start, amp_end, primers[0].pool)

    with open(output_amplicon_bed_path, "w") as output_amplicon_bed:
        for amplicon, coordinates in amplicons.items():
            amplicon_name = f"{schema_name}_{amplicon}"
            amplicon_bed_row = [
                reference_name,
                str(coordinates.start),
                str(coordinates.end),
                amplicon_name,
                str(coordinates.pool),
                "+"
            ]
            print("\t".join(amplicon_bed_row), file=output_amplicon_bed)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Convert a TSV file of primers into a FASTA file and two BED files (one for primers, one for amplicons).")
    parser.add_argument("input_tsv", type=argparse.FileType("r"), help="TSV file of primers.")
    parser.add_argument("--schema_name", type=str, help="Schema name (for FASTA output)")
    parser.add_argument("--reference", type=str, help="Reference name (for BED outputs)")

    args = parser.parse_args()
    schema_to_fasta(args.input_tsv, args.schema_name)
    schema_to_beds(args.input_tsv, args.schema_name, args.reference)
