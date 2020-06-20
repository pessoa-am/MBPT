#!/usr/bin/env python3

import sys
import xml.etree.ElementTree as ET
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an INSDSeq NCBI file to FASTA format.")
    parser.add_argument(
        "input_file", metavar="INSDXML_NCBI_FILE",
        type=str, help="Input INSDXML NCBI file (*.gbc)")
    parser.add_argument(
        "-cds", action="store_true",
        help='Retrieves the CDS portion of all cDNA sequences (if CDS '
             'information is available), adds "-" to correct the ORF, '
             'if necessary')
    parser.add_argument("-s", action="store_true",
                        help="Sorts sequences by name")
    parser.add_argument("-o", metavar="DESTINATION_FILE", type=str,
                        help="Output FASTA to file instead of stdout")
    return parser.parse_args()


def insd2seqrecords(input_file, cds=False, sort=False):
    locus = {}
    root = ET.parse(input_file).getroot()

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
            locus[locus_id]['cds'] = feature_data.findtext(
                "INSDFeature_location").strip()
            cs_data = feature_data.find(
                ".//INSDQualifier[INSDQualifier_name='codon_start']")
            if cs_data is not None:
                locus[locus_id]['cs'] = int(cs_data.findtext(
                    "INSDQualifier_value").strip())

        definition_data = seq_data.find(".//INSDSeq_definition")
        if definition_data is not None:
            locus[locus_id]['name'] = definition_data.text.strip().replace(
                locus_id, "").strip()

    if cds:
        for locus_id, locus_data in locus.items():
            if locus_data['cds'] is not None:
                start, end = locus_data['cds'].split("..")
                start, end = start.strip("<"), end.strip(">")
                start, end = int(start) - 1, int(end) - 1
                locus_data['seq'] = locus_data['seq'][start:end + 1]
                if locus_data['cs'] == 2:
                    locus_data['seq'] = f"--{locus_data['seq']}"
                elif locus_data['cs'] == 3:
                    locus_data['seq'] = f"-{locus_data['seq']}"

    if sort:
        sorted_loci = sorted(
            locus.keys(), key=lambda x: locus[x]['name'].upper())
    else:
        sorted_loci = locus.keys()

    seq_records = []
    for locus_id in sorted_loci:
        fasta_name = [locus_id, locus[locus_id]['name']]
        if cds and locus[locus_id]['cds'] is not None:
            fasta_name.append(locus[locus_id]['cds'])
        else:
            fasta_name.append(f"1..{locus[locus_id]['length']}")
        seq_records.append(SeqRecord(Seq(locus[locus_id]['seq']),
                           id="|".join(fasta_name), description=""))

    return seq_records


def main():
    args = parse_args()

    seq_records = insd2seqrecords(args.input_file, args.cds, args.s)

    if args.o:
        SeqIO.write(seq_records, args.o, "fasta")
    else:
        SeqIO.write(seq_records, sys.stdout, "fasta")


if __name__ == '__main__':
    main()
