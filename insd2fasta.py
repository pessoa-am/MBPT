#!/usr/bin/env python3

import sys
import xml.etree.ElementTree as ET

if len(sys.argv) < 2 or sys.argv[1] == "--help":
	print("Usage:\n\tinsd2fasta.py [INSDXML NCBI File (*.gbc)] [-cds] [-s] >[Destination FASTA File]\n"
		"Options:\n"
		"\t-cds\n"
		"\t\tRetrieves the CDS portion of all cDNA sequences (if CDS information is available), adds \"-\" to correct the ORF (if necessary and if ORF information is available) and is optional.\n"
		"\t-s\n"
		"\t\tSorts sequences by name and is optional.\n")
	sys.exit()

locus = {}
root = ET.parse(sys.argv[1]).getroot()

for seq_data in root.findall(".//INSDSeq"):
	locus_id = seq_data.findtext("INSDSeq_locus")
	locus[locus_id] = {
		'cds': None,
		'cs': None,
		'name': None,
		'seq': seq_data.findtext("INSDSeq_sequence").strip(),
		'length': seq_data.findtext("INSDSeq_length").strip()
	}

	feature_data = seq_data.find(".//INSDFeature[INSDFeature_key='CDS']")
	if feature_data is not None:
		locus[locus_id]['cds'] = feature_data.findtext("INSDFeature_location").strip()
		cs_data = feature_data.find(".//INSDQualifier[INSDQualifier_name='codon_start']")
		if cs_data is not None:
			locus[locus_id]['cs'] = int(cs_data.findtext("INSDQualifier_value").strip())

	definition_data = seq_data.find(".//INSDSeq_definition")
	if definition_data is not None:
		locus[locus_id]['name'] = definition_data.text.strip().replace(locus_id, "").strip()

if "-cds" in sys.argv:
	for locus_id, locus_data in locus.items():
		if locus_data['cds'] is not None:
			start, end = locus_data['cds'].split("..")
			start, end = start.strip("<"), end.strip(">")
			start, end = int(start) - 1, int(end) - 1
			locus_data['seq'] = locus_data['seq'][start:end+1]
			if locus_data['cs'] == 2:
				locus_data['seq'] = f"--{locus_data['seq']}"
			elif locus_data['cs'] == 3:
				locus_data['seq'] = f"-{locus_data['seq']}"

if "-s" in sys.argv:
	sorted_loci = sorted(locus.keys(), key=lambda x: locus[x]['name'].upper())
else:
	sorted_loci = locus.keys()

for locus_id in sorted_loci:
	fasta_name = [locus_id, locus[locus_id]['name']]
	if "-cds" in sys.argv and locus[locus_id]['cds'] is not None:
		fasta_name.append(locus[locus_id]['cds'])
	else:
		fasta_name.append(f"1..{locus[locus_id]['length']}")
	print(">" + "|".join(fasta_name))
	print(locus[locus_id]['seq'])


