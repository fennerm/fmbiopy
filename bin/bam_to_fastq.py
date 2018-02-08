#!/usr/bin/env python
"""Convert a bam file into forward, reverse and unpaired fastq files

Outputs: *.R1.fastq *.R2.fastq and *.unpaired.fastq

Unmapped reads are included in the outputs.

Usage:
  bam_to_fastq.py [--threads=N] BAM
  bam_to_fastq.py (-h | --help)

Options:
  -h --help     Show this screen
  --threads=N   Number of threads [default: 1]
"""
from docopt import docopt
from plumbum import local
from plumbum.cmd import (
    gzip,
    picard,
    samtools,
)


def extract_paired_reads(bam, threads='1'):
    tmp_bam = bam.with_suffix('.paired.bam')
    tmp_index = local.path(tmp_bam + ".bai")
    ((samtools['view', '-bh', '-@', threads, '-f', '1', bam] |
        samtools['sort']) > tmp_bam)()
    (samtools['index', tmp_bam])()
    fwd_fastq = bam.with_suffix(".R1.fastq")
    rev_fastq = bam.with_suffix(".R2.fastq")
    (picard['SamToFastq', 'I=' + tmp_bam, 'FASTQ=' + fwd_fastq,
            'SECOND_END_FASTQ=' + rev_fastq])()
    tmp_bam.delete()
    tmp_index.delete()
    return (fwd_fastq, rev_fastq)


def extract_unpaired_reads(bam, threads='1'):
    tmp_bam = bam.with_suffix('.unpaired.bam')
    tmp_index = local.path(tmp_bam + ".bai")
    (samtools['view', '-bh', '-@', threads, '-F', '1', bam] |
        samtools['sort'] > tmp_bam)()
    (samtools['index', tmp_bam])()
    unpaired_fastq = bam.with_suffix(".unpaired.fastq")
    (picard['SamToFastq', 'I=' + tmp_bam, 'FASTQ=' + unpaired_fastq])()
    tmp_bam.delete()
    tmp_index.delete()
    return unpaired_fastq


def main(bam, threads):
    (fwd_reads, rev_reads) = extract_paired_reads(bam, threads)
    unpaired_reads = extract_unpaired_reads(bam, threads)
    ps = [gzip.popen(fastq) for fastq in [fwd_reads, rev_reads, unpaired_reads]]
    for p in ps:
        p.wait()


if __name__ == '__main__':
    opt = docopt(__doc__)
    bam = local.path(opt['BAM'])
    threads = opt['--threads']
    main(bam, threads)
