#!/usr/bin/env python
"""Convert a sam file to a bam file"""
import sys

from plumbum import local
from plumbum.cmd import samtools

def sam_to_bam(sam, bam):
    """Convert a sam file to a bam file"""
    if sam.suffix != '.sam':
        raise OSError('File is not a .sam file')
    bam = sam.with_suffix('.bam')

    # Convert sam to sorted bam
    (samtools['view', '-Sb', str(sam)] | samtools['sort'] > str(bam))()

    # Index bam
    (samtools['index', str(bam)])()

if __name__ == '__main__':
    in_sam = local.path(sys.argv[1])
    if len(sys.argv) > 2:
        out_bam = local.path(sys.argv[2])
    else:
        out_bam = in_sam.with_suffix('.bam')
    sam_to_bam(in_sam, out_bam)
