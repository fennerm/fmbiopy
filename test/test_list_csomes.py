from plumbum import (
    local,
    ProcessExecutionError,
)
from pytest import raises

from fmbiopy.fmsystem import capture_stdout


def test_list_csomes_without_unmapped(trimmed_bam):
    list_csomes = local['bin/list_csomes']
    expected_output = ['simulated_fasta_500_' + str(num) for num in range(7)]
    assert capture_stdout(list_csomes[trimmed_bam]) == expected_output


def test_list_csomes_with_unmapped(bam_with_orphans):
    list_csomes = local['bin/list_csomes']
    expected_output = ['simulated_fasta_500_' + str(num) for num in range(7)]
    assert capture_stdout(list_csomes[bam_with_orphans]) == expected_output


def test_list_csomes_for_empty_bam(empty_bam):
    list_csomes = local['bin/list_csomes']
    with raises(ProcessExecutionError):
        assert capture_stdout(list_csomes[empty_bam]) == []
