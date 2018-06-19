#!/usr/bin/env python3
from collections import OrderedDict
import os
import sys
from tempfile import NamedTemporaryFile
from uuid import uuid4

import click
from pandas import DataFrame, read_csv
from plumbum import FG
from plumbum.cmd import bedtools

def endpoints_are_within_csome(bed_table, faidx):
    csome_length = faidx[faidx["CHROM"] == bed_table["CHROM"]]["LENGTH"]
    return bed_table["r1"] >= 0) and bed_table["r2"] <= csome_length.values[0]


@click.command()
@click.option("-r", "--reference", type=str, help="Reference sequence (fasta)")
@click.option(
    "-n", "--num-flanking", type=int, help="Number of flanking bases to extract"
)
@click.option(
    "--pre", type=int, help="Number of bases before variant to extract"
)
@click.option(
    "--post", type=int, help="Number of bases after variant to extract"
)
@click.argument("tsv", nargs=1)
def extract_flanking_bases(reference, pre, post, num_flanking, tsv):
    """Extract flanking sequence from a variant table."""
    if sum([bool(pre), bool(post), bool(num_flanking)]) > 1:
        sys.exit("Only one of -pre/-post/-n should be defined")
    variants = read_csv(tsv, sep="\t")
    if not len(variants):
        sys.stdout.write("\n")
        sys.exit()
    if num_flanking:
        r1 = variants["POS"] - num_flanking
        r2 = variants["POS"] + num_flanking
    elif pre:
        r1 = variants["POS"] - pre
        r2 = variants["POS"]
    elif post:
        r1 = variants["POS"]
        r2 = variants["POS"] + post
    bedfile_contents = DataFrame(
        OrderedDict([("CHROM", variants["CHROM"]), ("r1", r1), ("r2", r2)])
    )
    faidx = read_csv(
        reference + ".fai",
        names=["CHROM", "LENGTH"],
        sep="\t",
        usecols=[0, 1])
    has_valid_endpoints = bedfile_contents.apply(
        endpoints_are_within_csome,
        axis=1,
        faidx=faidx
    )
    bedfile_contents = bedfile_contents[has_valid_endpoints]
    tmpfile = os.path.join("/tmp", uuid4().hex, "flanking.bed")
    with NamedTemporaryFile() as tmpfile:
        bedfile_contents.to_csv(
            tmpfile.name, sep="\t", header=False, index=False
        )
        bedtools["getfasta", "-fi", reference, "-bed", tmpfile.name] & FG


if __name__ == "__main__":
    extract_flanking_bases()
