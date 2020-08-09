# MBPT (Molecular Biology Python Tools)

This is a collection of python scripts I've been making over the years to assist me in molecular biology and bioinformatics related work. Some were originally written in Perl, but I've since ported them to Python3. GPL-3 applies unless otherwise stated.

## sortfasta.py

Sorts the ids of a FASTA file by name.

Requires `argparse` and `biopython`.

## insd2fasta.py

Converts an INSDSeq NCBI file to FASTA format. Has the possibility of retrieving the CDS portion of all cDNA sequences (if CDS information is available), adds \"-\" to correct the ORF if necessary. It also sorts the ids of a FASTA file by name.

Requires `argparse` and `xml.etree.ElementTree`.

## blastxml2fasta.py

A simple BLAST XML to FASTA converter. Adds 5' gaps to facilitate downstream alignment with query sequence.

Requires `argparse` and `biopython`.

## consensus.py

Generates a consensus sequence from a pre-aligned FASTA file.

Requires `argparse` and `biopython`.

## pdw2fasta.py

Converts a pDRAW32 file into FASTA.

Requires `argparse` and `biopython`.

## pdw2json.py

Converts a pDRAW32 file into JSON.

Requires `argparse`.

## get_misprime_libs.py

Converts mispriming libraries downloaded directly from the Primer3Plus repository into Python a file that can be imported accordingly. Particularly useful for [primer3-py](https://github.com/libnano/primer3-py). By default creates a file called `misprime_libs.py`.

Requires `argparse` and `biopython`.
