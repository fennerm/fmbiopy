#!/usr/bin/env python
"""Align reads with bowtie2 and sort
print(sys.argv)

This is a wrapper around bowtie2 which handles the indexing step before
alignment, and converts the output sam to a sorted and indexed bam file.

Output is sent to stdout. -S argument should not be passed since the stdout
of Bowtie2 needs to be piped to samtools.

Requires: bowtie2, samtools

Usage:
  align_and_sort.py [args]*
  align_and_sort.py (-h | --help)

Options:
  -h --help                 Show this screen
"""
from __future__ import print_function
import sys
from uuid import uuid4

from plumbum import (
    BG,
    local,
)
from plumbum.cmd import (
    bowtie2,
    samtools,
)


def _get_reference(args):
    """Get the path to the reference seq from a list of bowtie2 arguments"""
    ref_index = _optget(args, "-x")

    reference = local.path(ref_index + '.fa')
    if not reference.exists():
        reference = local.path(ref_index + '.fasta')
        if not reference.exists():
            raise ValueError('Fasta not found')
    return reference


def _get_threads(args):
    """Get the number of threads from a list of bowtie2 arguments"""
    try:
        threads = _optget(args, "-p", "--threads")
    except ValueError:
        threads = 1
    return threads

def _optget(args, short, long=None):
    """Parse command line arguments to fetch the value of an option

    Parameters
    ----------
    args: List[str]
        Command line arguments
    short: str
        The short name of the option (including dash)
    long: str, optional
        The long name of the option (including dashes)

    Returns
    -------
    str
        The option value.

    Raises
    ------
    ValueError
        If option is not present in args
    """
    try:
        # Try parse as a short option
        opt = args[args.index(short) + 1]
    except ValueError:
        if long:
            # Try parse as a long option
            for arg in args:
                if arg.startswith(long):
                    if arg == long:
                        opt = args[args.index(arg) + 1]
                    else:
                        opt = arg.split("=")[1]
            try:
                opt
            except NameError:
                raise ValueError("Opt not present in args")
        raise
    return opt

def align_and_sort(args):
    ref = _get_reference(args)
    fai_indices = local.path(ref + '.fai')
    if not fai_indices:
        print('Indexing fasta with samtools', sys.stderr)
        samtools_index = samtools('faidx', ref) & BG
        samtools_index.wait()

    bowtie2_index = ref.with_suffix('.1.bt2')
    if not bowtie2_index.exists():
        print('Indexing fasta with bowtie2')
        prefix = ref.with_suffix('')
        build = local['bowtie2-build']
        p = build[ref, prefix] & BG
        p.wait()

    threads = _get_threads(args)
    align_pipe = (
        bowtie2.__getitem__(args) |
        samtools['view', '-bhS'] |
        samtools['sort', '-@', threads, '-T', uuid4().hex[1:8]])
    p = align_pipe & BG(stdout=sys.stdout, stderr=sys.stderr)
    p.wait()


if __name__ == '__main__':
    if sys.argv[1] in ['-h', '--help']:
        print(__doc__)
        sys.exit(0)
    align_and_sort(sys.argv[1:])
