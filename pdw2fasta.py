#!/usr/bin/env python3
import sys
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args():
	parser = argparse.ArgumentParser(description="Converts a pDRAW32 file into FASTA")
	parser.add_argument("input_file", metavar="PDW_FILE", help="Input pDRAW32 file (*.pdw)")
	parser.add_argument("-o", metavar="OUTPUT_FILE", help="Output FASTA to file instead of stdout")

	return parser.parse_args()

def pdw2seqrecord(input_file):
	
	with open(input_file, "r") as f:
		content = f.readlines()

	content = [x.strip().split(maxsplit=1) for x in content]

	seq = ""
	name = ""
	
	for line in content:
		if line[0].isdigit():
			seq += line[1].replace(" ", "")
		if line[0] == "DNAname":
			name = line[1]

	output = SeqRecord(Seq(seq), name, "", "")
	return output


def main():
	args = parse_args()

	seq_record = pdw2seqrecord(args.input_file)
	
	if args.o:
		SeqIO.write(seq_record, args.o, "fasta")
	else:
		SeqIO.write(seq_record, sys.stdout, "fasta")

if __name__ == "__main__":
	main()
