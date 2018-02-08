from fmbiopy.fmtest import assert_script_produces_files
from plumbum.cmd import samtools


def test_bam_to_fastq(indexed_small_bam):
    suffixes = [".R1.fastq.gz", ".R2.fastq.gz", ".unpaired.fastq.gz"]
    expected_output = [indexed_small_bam.with_suffix(suff) for suff in suffixes]
    expected_outdir = indexed_small_bam.dirname

    assert_script_produces_files("bin/bam_to_fastq.py", indexed_small_bam,
                                 expected_output, outdir=expected_outdir)
