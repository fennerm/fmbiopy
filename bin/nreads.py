#!/usr/bin/env python
"""Calculate the number of reads in a fastq or fastq.gz file

Usage:
    nreads.py FILE
    nreads.py (-h | --help)
"""
from __future__ import print_function

from plumbum import local

from fmbiopy.fmbio import count_reads
from fmbiopy.fmparse import helpful_docopt

if __name__ == "__main__":
    opts = helpful_docopt(__doc__)
    filename = local.path(opts['FILE'])
    print(count_reads(filename))
