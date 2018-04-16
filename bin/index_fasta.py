#!/usr/bin/env python
"""Index a .fasta file using samtools, bowtie2 and picard

Usage:
  index_fasta.py FASTA
  index_fasta.py (-h --help)
"""
from docopt import docopt
from plumbum import (
    FG,
    local,
    )
from plumbum.cmd import (
        picard,
        samtools,
        )

def main(fasta):
    fasta = local.path(fasta)

    fai = local.path(fasta + '.fai')
    if not fai.exists():
        samtools['faidx', fasta] & FG

    pic_idx = fasta.with_suffix('.dict')
    if not pic_idx.exists():
        picard['CreateSequenceDictionary', 'R=' + fasta, 'O=' + pic_idx] & FG

    bt2_idx = fasta.with_suffix('.1.bt2')
    if not bt2_idx.exists():
        bowtie2_build = local['bowtie2-build']
        bowtie2_build[fasta, fasta.with_suffix('')] & FG

if __name__ == "__main__":
    opts = docopt(__doc__)
    main(opts['FASTA'])
