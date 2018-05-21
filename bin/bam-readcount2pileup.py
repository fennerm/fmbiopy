#!/usr/bin/env python
"""Convert output from bam-readcount to a allele count table.

The primary benefit of this transformation is that the output will have an even
number of columns and so will be easier to parse for later operations.

Much of the information from the readcount file is discarded, leaving only raw
counts of each allele.
"""
import csv
import sys

from plumbum import local


def count_columns(tsv_file):
    return int(local["max_ncol.awk"](tsv_file).rstrip())


def get_allele_counts(allele_pileup):
    """Given an allele pileup from the readcount file, extract the allele name
    and count as a string.

    E.g get_allele_counts('T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0') = 'T:0'
    """
    return ":".join(allele_pileup.split(":")[0:2])


def parse_row(row, num_columns):
    """Convert row from input tsv to output tsv format."""
    counts = [get_allele_counts(pile) for pile in row[5:]]
    num_missing_columns = num_columns - 3 - len(counts)
    return [row[0], row[1], row[3]] + counts + (["NA"] * num_missing_columns)


def print_header(num_columns):
    """Print the output file header to stdout."""
    count_names = ["count" + str(i) for i in range(num_columns - 3)]
    print("id\tpos\tcov\t" + "\t".join(count_names))


def bamreadcount2pileup(input_tsv):
    num_columns = count_columns(input_tsv)
    print_header(num_columns)
    with open(input_tsv, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            print("\t".join(parse_row(row, num_columns)))


if __name__ == "__main__":
    bamreadcount2pileup(sys.argv[1])
