#!/usr/bin/env python3

import argparse
import urllib.request
from io import TextIOWrapper
from Bio import SeqIO

# Copyright 2020 Alberto Pessoa <pessoa.am@gmail.com>
#
# See COPYING for full GPLv3 license.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# # You should have received a copy of the
# GNU General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.

lib_license = """# Copyright (c) 2006 - 2018
# by Andreas Untergasser and Harm Nijveen
# All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# # You should have received a copy of the
# GNU General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.
"""


def main():
    parser = argparse.ArgumentParser(
        description='Converts mispriming libraries downloaded directly from '
                    'the Primer3Plus repository into Python a file that can '
                    'be imported accordingly')
    parser.add_argument(
        '-o', '--output-file', default='misprime_libs.py',
        help='The name of the Python file to store the output dictionaries '
             '(default: misprime_libs.py)')
    args = parser.parse_args()

    url_dict = {'HUMAN': 'https://raw.githubusercontent.com/primer3-org/primer3plus/main/server/mispriming_lib/humrep_and_simple.txt',
                'RODENT_AND_SIMPLE': 'https://raw.githubusercontent.com/primer3-org/primer3plus/main/server/mispriming_lib/rodrep_and_simple.txt',
                'RODENT': 'https://raw.githubusercontent.com/primer3-org/primer3plus/main/server/mispriming_lib/rodent_ref.txt',
                'DROSOPHILA': 'https://raw.githubusercontent.com/primer3-org/primer3plus/main/server/mispriming_lib/drosophila_w_transposons.txt'
                }

    misprime_libs = {}

    for lib, url in url_dict.items():
        with urllib.request.urlopen(url) as response:
            fasta = SeqIO.parse(
                TextIOWrapper(response, encoding='utf-8'), 'fasta')
            misprime_libs[lib] = {record.description:
                                  str(record.seq) for record in fasta}

    with open(args.output_file, 'w') as pyfile:
        pyfile.write(f'{lib_license}\n\n')
        for lib, lib_dict in misprime_libs.items():
            pyfile.write(f'{lib} = ' + '{')
            first = True
            content = []
            for key, value in lib_dict.items():
                key = key.replace("\\", "_")
                tab = ''
                if not first:
                    tab = ' ' * len(lib) + '    '
                first = False
                content.append(f'{tab}"{key}": "{value}"')
            pyfile.write(',\n'.join(content) + '\n')
            pyfile.write(' ' * len(lib) + '    }\n')


if __name__ == '__main__':
    main()
