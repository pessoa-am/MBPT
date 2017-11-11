#!/usr/bin/env python3

import sys
import xml.etree.ElementTree as ET

if len(sys.argv) < 2 or sys.argv[1] == "--help":
	sys.exit("Usage:\nblastxml2fasta.py [BLAST XML File (*.xml)] >[Destination FASTA File]")

hits = {}
root = ET.parse(sys.argv[1]).getroot()

for hit in root.findall(".//Hit"):
	hit_def = hit.find("Hit_def").text.strip()
	hit_accession = hit.find("Hit_accession").text.strip()
	locus = f"{hit_accession}|{hit_def}"
	hsp = hit.find("Hit_hsps/Hsp")
	cs = int(hsp.find("Hsp_query-from").text)
	seq = hsp.find("Hsp_hseq").text.strip()

	for i in range(2, cs + 1):
		seq = "-" + seq

	hits[locus] = seq

for locus in hits:
	print(f">{locus}")
	print(hits[locus])
