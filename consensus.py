#!/usr/bin/env python3
import sys
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a consensus sequence from a FASTA file")
    parser.add_argument(
        "input_file", metavar="FASTA_FILE", help="Input FASTA file")
    parser.add_argument(
        "-o", metavar="OUTPUT_FILE",
        help="Output FASTA to file instead of stdout")

    return parser.parse_args()


def consensus(seq_records, consensus_id='Consensus'):

    matrix = [
        ["-", "a", "c", "g", "t", "r", "y", "s", "w", "k", "m", "b", "d", "h", "v", "n"],
        ["a", "A", "M", "R", "W", "R", "H", "V", "W", "D", "M", "N", "D", "H", "V", "N"],
        ["c", "M", "C", "S", "Y", "V", "Y", "S", "H", "B", "M", "B", "N", "H", "V", "N"],
        ["g", "R", "S", "G", "K", "R", "B", "S", "D", "K", "V", "B", "D", "N", "V", "N"],
        ["t", "W", "Y", "K", "T", "D", "Y", "B", "W", "K", "H", "B", "D", "H", "N", "N"],
        ["r", "R", "V", "R", "D", "R", "N", "V", "D", "D", "V", "N", "D", "N", "V", "N"],
        ["y", "H", "Y", "B", "Y", "N", "Y", "B", "H", "B", "H", "B", "N", "H", "N", "N"],
        ["s", "V", "S", "S", "B", "V", "B", "S", "N", "B", "V", "B", "N", "N", "V", "N"],
        ["w", "W", "H", "D", "W", "D", "H", "N", "W", "D", "H", "N", "D", "H", "N", "N"],
        ["k", "D", "B", "K", "K", "D", "B", "B", "D", "K", "N", "B", "D", "N", "N", "N"],
        ["m", "M", "M", "V", "H", "V", "H", "V", "H", "N", "M", "N", "N", "H", "V", "N"],
        ["b", "N", "B", "B", "B", "N", "B", "B", "N", "B", "N", "B", "N", "N", "N", "N"],
        ["d", "D", "N", "D", "D", "D", "N", "N", "D", "D", "N", "N", "D", "N", "N", "N"],
        ["h", "H", "H", "N", "H", "N", "H", "N", "H", "N", "H", "N", "N", "H", "N", "N"],
        ["v", "V", "V", "V", "N", "V", "N", "V", "N", "N", "V", "N", "N", "N", "V", "N"],
        ["n", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N"]
    ]

    target = []
    for record in seq_records:
        target.append(str(record.seq))

    le = min([len(i)-1 for i in target])

    contarget = ""
    j = 0

    for i in range(le + 1):
        for sequence in target:
            nuc = sequence[i]
            if j == 0:
                con = nuc
            if re.match("[acgtryswkmhdvnb]", con) and j != 0:
                case = 1
            else:
                case = 0
            nuc = nuc.translate(str.maketrans(
                '-ACGTRYSWKMBDHVNacgtryswkmbdhvn',
                '0123456789ABCDEF123456789ABCDEF'))
            con = con.translate(str.maketrans(
                '-ACGTRYSWKMBDHVNacgtryswkmbdhvn',
                '0123456789ABCDEF123456789ABCDEF'))
            nuc = int(nuc, 16)
            con = int(con, 16)
            con = matrix[nuc][con]
            if case == 1:
                con = con.lower()
            j += 1
        j = 0
        contarget += con

    output = SeqRecord(Seq(contarget), consensus_id, "", "")
    return output


def main():
    args = parse_args()

    seq_records = []

    for record in SeqIO.parse(args.input_file, "fasta"):
        seq_records.append(record)

    record = consensus(seq_records)

    if args.o:
        SeqIO.write(record, args.o, "fasta")
    else:
        SeqIO.write(record, sys.stdout, "fasta")


if __name__ == "__main__":
    main()
