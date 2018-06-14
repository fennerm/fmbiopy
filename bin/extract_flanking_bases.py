#!/usr/bin/env python
import os
from tempfile import NamedTemporaryFile
from uuid import uuid4

import click
from pandas import DataFrame, read_csv
from plumbum import FG
from plumbum.cmd import bedtools


@click.command()
@click.option("-r", "--reference", help="Reference sequence (fasta)")
@click.option(
    "-n", "--num-flanking", type=int, help="Number of flanking bases to extract"
)
@click.argument("tsv", nargs=1)
def extract_flanking_bases(reference, num_flanking, tsv):
    """Extract flanking sequence from a variant table."""
    variants = read_csv(tsv, sep="\t")
    bedfile_contents = DataFrame(
        {
            "CHROM": variants["CHROM"],
            "r1": variants["POS"] - num_flanking,
            "r2": variants["POS"] + num_flanking,
        }
    )
    bedfile_contents = bedfile_contents[["CHROM", "r1", "r2"]]
    tmpfile = os.path.join("/tmp", uuid4().hex, "flanking.bed")
    with NamedTemporaryFile() as tmpfile:
        bedfile_contents.to_csv(
            tmpfile.name, sep="\t", header=False, index=False
        )
        bedtools["getfasta", "-fi", reference, "-bed", tmpfile.name] & FG


if __name__ == "__main__":
    extract_flanking_bases()
