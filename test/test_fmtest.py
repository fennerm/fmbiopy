"""Test suite for fmbiopy.fmtest"""
from uuid import uuid4

from plumbum import local
from plumbum.cmd import samtools

from fmbiopy.fmbio import count_reads
from fmbiopy.fmsystem import capture_stdout


def test_cd(cd, startdir):
    assert local.cwd != startdir


def test_bam_with_orphans(sandbox, bam_with_orphans, trimmed_bam):
    fastqs = [sandbox / (uuid4().hex + suffix)
              for suffix in ['.R1.fastq', '.R2.fastq', '.unpaired.fastq']]
    # Check that naively converting to fastq leads to uneven number of pairs
    samtools['fastq', '-1', fastqs[0], '-2', fastqs[1], '-0', fastqs[2],
             bam_with_orphans]()
    assert count_reads(fastqs[0]) != count_reads(fastqs[1])
