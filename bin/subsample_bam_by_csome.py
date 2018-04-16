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

from docopt import docopt
from plumbum.cmd import (
    extract_csome,
    list_csomes,
    samtools,
)

from fmbiopy.fmlist import exclude_blank


def main(n, inbam, outbam):
    n = int(n)
    csome_names = exclude_blank(list_csomes(inbam).split("\n"))
    random_csomes = sample(set(csome_names), n)
    extract_args = [inbam] + random_csomes
    (extract_csome.__getitem__(extract_args) > outbam)()
    sys.exit()


if __name__ == "__main__":
    opts = docopt(__doc__)
    main(n=opts["--ncsomes"], inbam=opts["INBAM"], outbam=opts["--output"])
