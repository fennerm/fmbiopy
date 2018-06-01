#!/usr/bin/env python
"""Calculate the number of reads in a fastq or fastq.gz file

Usage:
    nreads.py FILE
    nreads.py (-h | --help)
"""
from __future__ import print_function

from docopt import docopt
from plumbum import local

from fmbiopy.bio import count_reads

if __name__ == "__main__":
    opts = docopt(__doc__)
    filename = local.path(opts['FILE'])
    print(count_reads(filename))
