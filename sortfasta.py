#!/usr/bin/env python3
import sys

if len(sys.argv) < 2 or sys.argv[1] == "--help":
	sys.exit("Usage:\nsortfasta.py [FASTA Sequence File] >[Destination FASTA File]\n")

target = {}
name = ""
with open(sys.argv[1], "r") as FASFILE:
	for line in FASFILE:
		line = line.strip()
		if line.startswith(">"):
			name = line.replace(">", "")
			target[name] = ""
		else:
			target[name] += line

sorted_names = sorted(target.keys(), key=lambda x: x.upper())

for name in sorted_names:
	print(f">{name}")
	print(target[name])
