#!/usr/bin/env python3
import sys
import argparse
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sort the ids of a FASTA file by name")
    parser.add_argument(
        "input_file", metavar="FASTA_FILE", help="Input FASTA file")
    parser.add_argument(
        "-o", metavar="OUTPUT_FILE",
        help="Output FASTA to file instead of stdout")

    return parser.parse_args()


def sort_seqrecords(seq_records):
    records_dict = {}

    for record in seq_records:
        records_dict[record.id] = record

    sorted_ids = sorted(records_dict.keys(), key=lambda x: x.upper())

    seq_records = []

    for key in sorted_ids:
        seq_records.append(records_dict[key])

    return seq_records


def main():
    args = parse_args()

    seq_records = []

    for record in SeqIO.parse(args.input_file, "fasta"):
        seq_records.append(record)

    seq_records = sort_seqrecords(seq_records)

    if args.o:
        SeqIO.write(seq_records, args.o, "fasta")
    else:
        SeqIO.write(seq_records, sys.stdout, "fasta")


if __name__ == "__main__":
    main()
