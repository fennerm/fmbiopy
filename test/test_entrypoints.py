"""Test module CLI entrypoints."""
from plumbum import local

from fmbiopy.system import capture_stdout


def test_nreads(fasta):
    assert capture_stdout(local["bin/nseqs"][fasta["fasta"]]) == ["50"]
