#!/usr/bin/env python
"""Calculate the number of reads in a fastq or fastq.gz file

Usage:
    nreads.py FILE
    nreads.py (-h | --help)
"""
from __future__ import print_function

from fmbiopy.fmparse import helpful_docopt
from plumbum import local
from plumbum.cmd import (
    bc,
    echo,
    wc,
    zcat,
)


def nreads(path):
    """Calculate the number of reads in a fastq or fastq.gz file

    Results written to stdout

    Parameters
    ----------
    path: pathlike
        A fastq or fastq.gz file

    Returns
    -------
    int
        The number of reads in the file
    """
    path = local.path(path)
    if path.suffix == '.gz':
        nreads = zcat[path] | wc['-l']
        n = int(int((nreads.run()[1]).split(' ')[0]) / 4)
    else:
        n = int(int((wc.run(['-l', path])[1]).split(' ')[0]) / 4)

    return n


if __name__ == "__main__":
    opts = helpful_docopt(__doc__)
    print(nreads(opts['FILE']))
