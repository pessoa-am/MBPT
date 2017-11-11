#!/usr/bin/env python3

import sys
import re

if len(sys.argv) < 2 or sys.argv[1] == "--help":
	sys.exit("Usage:\npdw2fasta.py [pDRAW32 File (*.pdw)] >[Destination FASTA File]")

with open(sys.argv[1], "r") as f:
	content = f.readlines()

content = [x.strip().split(maxsplit=1) for x in content]

seq = ""

for line in content:
	if line[0].isdigit():
		seq += line[1].replace(" ", "")
	if line[0] == "DNAname":
		print(f">{line[1]}")

print(seq)
