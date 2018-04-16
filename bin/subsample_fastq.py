#!/usr/bin/env python
"""Subsample reads from a pair of fastq or fastq.gz files using reformat.sh

Usage:
  subsample_fastq.py (--sample=<N>) IN1 IN2 OUT1 OUT2
  subsample_fastq.py (-h --help)

Inputs:
  IN1, IN2      Input fastq(.gz) files
  OUT1, OUT2    Output subsampled files

Options:
  -h --help         Show this screen
  -n --sample=<N>   Number of reads to sample
"""

from docopt import docopt
from plumbum import local


def subsample_fastq(in1, in2, out1, out2, n):
    """Subsample reads from a pair of fastq files using seqtk

    Parameters
    ----------
    in1, in2
        The input pair
    out1, out2
        The paths to the subsampled output files
    """
    reformat = local['reformat.sh']
    reformat[
            'in=' + in1, 'in2=' + in2, 'out=' + out1, 'out2=' + out2,
            'samplereadstarget=' + n]()


if __name__ == "__main__":
    opts = docopt(__doc__)

    subsample_fastq(opts['IN1'], opts['IN2'], opts['OUT1'], opts['OUT2'],
                    n=opts['--sample'])
