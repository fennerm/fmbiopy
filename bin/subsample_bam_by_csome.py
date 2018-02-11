#!/usr/bin/env python
"""Subsample N chromosomes/contigs from a bam file

Only chromosomes/contigs which have at least one alignment are sampled

Usage:
  subsample_bam_by_csome.py -n N -o OUTBAM INBAM

Options:
  -n, --ncsomes=N       Number of contigs to sample
  -o, --output=OUTBAM   Output file
"""
import sys
from random import sample

from plumbum.cmd import (
    extract_csome,
    list_nonempty_csomes,
    samtools,
)

from fmbiopy.fmparse import helpful_docopt
from fmbiopy.fmlist import exclude_blank


def main(n, inbam, outbam):
    n = int(n)
    csome_names = exclude_blank(list_nonempty_csomes(inbam).split("\n"))
    random_csomes = sample(set(csome_names), n)
    args = [inbam] + random_csomes
    (extract_csome.__getitem__(args) | samtools['sort'] > outbam)()


if __name__ == "__main__":
    opts = helpful_docopt(__doc__)
    main(n=opts["--ncsomes"], inbam=opts["INBAM"], outbam=opts["--output"])
