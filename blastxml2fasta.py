#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert BLAST XML to FASTA format. Adds 5' gaps to "
                    "facilitate downstream alignment with query sequence.")
    parser.add_argument(
        "input_file", metavar="XML_FILE", help="Input BLAST XML file (*.xml)")
    parser.add_argument(
        "-o", metavar="OUTPUT_FILE",
        help="Output FASTA to file instead of stdout")

    return parser.parse_args()


def blastxml2seqrecords(input_file):
    hits = {}
    blast_record = NCBIXML.read(open(input_file, 'r'))
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            hit_def = alignment.hit_def
            hit_accession = alignment.accession
            locus = f"{hit_accession}|{hit_def}"
            cs = int(hsp.query_start)
            seq = hsp.sbjct
            for i in range(2, cs + 1):
                seq = "-" + seq
            hits[locus] = seq
            break

    seq_records = []
    for locus in hits:
        seq_records.append(
            SeqRecord(Seq(hits[locus]), id=locus, description=""))
    return seq_records


def main():
    args = parse_args()

    seq_records = blastxml2seqrecords(args.input_file)

    if args.o:
        SeqIO.write(seq_records, args.o, "fasta")
    else:
        SeqIO.write(seq_records, sys.stdout, "fasta")
        pass


if __name__ == "__main__":
    main()
