#!/usr/bin/env python3

import sys
import re
import json
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Converts a pDRAW32 file into JSON")
    parser.add_argument(
        "input_file", metavar="PDW_FILE", help="Input pDRAW32 file (*.pdw)")
    parser.add_argument(
        "-o", metavar="OUTPUT_FILE",
        help="Output JSON to file instead of stdout")

    return parser.parse_args()


if len(sys.argv) < 2 or sys.argv[1] == "--help":
    sys.exit(
        "Usage:\npdw2json.py [pDRAW32 File (*.pdw)] >[Destination JSON File]")


def pdw2json(input_file):

    labels = ['VERSION', 'DNAname', 'IScircular', 'PlotLinear', 'Element',
              'AvailAll', 'Commercial', 'InStock', 'Fiveprime', 'Blunt',
              'Threeprime', 'Palindromic', 'NonPalindromic', 'Degenerated',
              'Nondegenerated',
              'Interrupted', 'Noninterrupted', 'Length', 'InEntire', 'MinCuts',
              'MaxCuts', 'InRegion', 'RegionStart', 'MinCutsRegion',
              'MaxCutsRegion', 'RegionEnd', 'ForceOnly', 'ForcedEnzyme',
              'HEADER', 'Sequence']

    with open(sys.argv[1], 'r') as f:
        content = f.readlines()

    content = [x.strip().split(maxsplit=1) for x in content]

    data = dict((x, '') for x in labels)

    for key in ['Element', 'InRegion', 'RegionStart', 'MinCutsRegion',
                'MaxCutsRegion', 'RegionEnd', 'ForcedEnzyme', 'HEADER']:
        data[key] = []

    for line in content:
        m = re.search(
            "(InRegion|RegionStart|MinCutsRegion|MaxCutsRegion|RegionEnd)\d?",
            line[0])
        if m:
            data[m.group(1)].append(int(line[1]))
        if line[0].isdigit():
            data['Sequence'] += line[1].replace(" ", "")
        if line[0] in labels:
            if line[0] not in ['VERSION', 'DNAname', 'Element', 'InRegion',
                               'RegionStart', 'MinCutsRegion', 'MaxCutsRegion',
                               'RegionEnd', 'ForcedEnzyme', 'HEADER',
                               'Sequence']:
                if line[1] == 'YES':
                    data[line[0]] = True
                elif line[1] == 'NO':
                    data[line[0]] = False
                else:
                    data[line[0]] = int(line[1])
            elif line[0] in ['VERSION', 'DNAname']:
                data[line[0]] = line[1]
            elif line[0] == 'Element':
                elements = line[1].split()
                i = len(elements) - 4
                for e in range(i, len(elements)):
                    elements[e] = int(elements[e])
                data[line[0]].append(elements)
            elif line[0] == 'ForcedEnzyme':
                data[line[0]].append(line[1])
            elif line[0] == 'HEADER':
                data[line[0]].append(line[1])

    return json.dumps(data, separators=(',', ':'))


def main():
    args = parse_args()

    output = pdw2json(args.input_file)

    if args.o:
        with open(args.o, "w") as JSONFILE:
            JSONFILE.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
