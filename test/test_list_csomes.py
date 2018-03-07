from plumbum import (
    local,
    ProcessExecutionError,
)
from pytest import (
    fixture,
    raises,
)

from fmbiopy.fmsystem import capture_stdout


@fixture
def expected_seqids():
    return ['simulated_fasta_500_' + str(num) for num in range(7)]


@fixture
def list_csomes():
    return local['bin/list_csomes']


def test_list_csomes_with_fasta_input(list_csomes, partial_fasta,
                                      expected_seqids):
    seqids = capture_stdout(list_csomes[partial_fasta['fasta']])
    assert seqids == expected_seqids


def test_list_csomes_without_unmapped(list_csomes, trimmed_bam,
                                      expected_seqids):
    assert capture_stdout(list_csomes[trimmed_bam]) == expected_seqids


def test_list_csomes_with_unmapped(list_csomes, bam_with_orphans,
                                   expected_seqids):
    list_csomes = local['bin/list_csomes']
    assert capture_stdout(list_csomes[bam_with_orphans]) == expected_seqids


def test_list_csomes_for_empty_bam(list_csomes, empty_bam):
    with raises(ProcessExecutionError):
        assert capture_stdout(list_csomes[empty_bam]) == []
