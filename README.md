# MBPT (Molecular Biology Python Tools)

This is a collection of python scripts I've been making over the years to assist me in molecular biology and bioinformatics related work. Some were originally written in Perl, but I've since ported them to Python3. GPL-3 applies unless otherwise stated.

## sortfasta.py

Sorts the ids of a FASTA file by name.

## insd2fasta.py

Converts an INSDSeq NCBI file to FASTA format. Has the possibility of retrieving the CDS portion of all cDNA sequences (if CDS information is available), adds \"-\" to correct the ORF if necessary. It also sorts the ids of a FASTA file by name.

Requires `xml.etree.ElementTree`.

## blastxml2fasta.py

A simple BLAST XML to FASTA converter. Adds 5' gaps to facilitate downstream alignment with query sequence.

Requires `xml.etree.ElementTree`.

## consensus.py

Generates a consensus sequence from a pre-aligned FASTA file.

## pdw2fasta.py

Converts a pDRAW32 file into FASTA.

## pdw2json.py

Converts a pDRAW32 file into JSON.

