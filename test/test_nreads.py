"""Test module for nreads script"""

from bin.nreads import nreads

def test_nreads(dat):
    input_fastq = dat['tiny']['fwd_reads'][0]
    gzipped_fastq = dat['tiny']['zipped_fwd_reads'][0]
    for f in [input_fastq, gzipped_fastq]:
        assert nreads(f) == 40
