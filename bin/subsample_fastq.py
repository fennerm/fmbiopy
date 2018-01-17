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

from biofile import Fastq
from plumbum import (
        local,
        LocalPath,
        )

from fmbiopy.fmparse import helpful_docopt
from fmbiopy.fmpaths import as_paths

def subsample_fastq(
        in1: Fastq,
        in2: Fastq,
        out1: LocalPath,
        out2: LocalPath,
        n: int)-> None:
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
    opts = helpful_docopt(__doc__)
    readkeys = ['IN1', 'IN2']
    inps = [Fastq(opts[key], determine_gzip=True) for key in readkeys]
    subopts = [opts['OUT1'], opts['OUT2']]
    outs = as_paths(subopts)

    subsample_fastq(inps[0], inps[1], outs[0], outs[1], n=opts['--sample'])
