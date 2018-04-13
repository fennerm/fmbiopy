#!/usr/bin/env python
"""Extract coverage info from a .vcf file and export to a text file."""
import os

import click


@click.command()
@click.option('--output', '-o', help='Output text file')
@click.argument('input')
def extract_variant_coverage(input, output):
    """Extract coverage info from a .vcf file and export to a text file."""
    os.system('egrep -v "^#" {} | \\'
              'cut -f 8 | \\'
              'sed \'s/^.*;DP=\\([0-9]*\\);.*$/\\1/\' > {}'.format(input,
                                                                   output))

if __name__ == '__main__':
    extract_variant_coverage()
