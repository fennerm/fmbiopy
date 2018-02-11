from uuid import uuid4

from Bio import SeqIO
from pytest import mark

from fmbiopy.fmbio import *
from fmbiopy.fmlist import none
from fmbiopy.fmpaths import (
    all_exist,
    apply_is_empty,
    is_empty,
)


def test_count_reads_output_correct_for_sample_bam(dat):
    input_bam = dat['tiny']['bam'][0]
    assert count_reads(input_bam) == 20


def test_count_reads_output_correct_for_sample_sam(dat):
    input_sam = dat['tiny']['sam'][0]
    assert count_reads(input_sam) == 16


def test_count_reads_output_correct_for_empty_file(gen_tmp):
    tmp = gen_tmp(empty=True, suffix='.fastq')
    assert count_reads(tmp) == 0


def test_merge_bams(dat, sandbox):
    bams = dat['tiny']['bam'][0:2]
    output_file = sandbox / (uuid4().hex + '.bam')
    merge_bams(bams, output_file)
    assert count_reads(output_file) == count_reads(
        bams[0]) + count_reads(bams[1])


def test_to_fastq(dat, sandbox):
    bam = dat['small']['bam'][0]
    output_prefix = sandbox / uuid4().hex
    suffixes = ['.R1.fastq.gz', '.R2.fastq.gz', '.unpaired.fastq.gz']
    expected_output = [local.path(output_prefix + suff) for suff in suffixes]
    to_fastq(bam, output_prefix)
    assert all_exist(expected_output) and none(apply_is_empty(expected_output))


@mark.parametrize("contig_length", [0, 10])
@mark.parametrize("num_contigs", [0, 2])
def test_simulate_fasta(sandbox, contig_length, num_contigs):
    simulated_fasta = sandbox / (uuid4().hex + '.fasta')
    simulate_fasta(num_contigs, contig_length, simulated_fasta)
    if num_contigs == 0:
        assert is_empty(simulated_fasta)
    else:
        for record in SeqIO.parse(simulated_fasta, "fasta"):
            assert record.id
            assert len(record.seq) == contig_length


@mark.parametrize("include_n", [True, False])
def test_simulate_fasta_includes_Ns(sandbox, include_n):
    simulated_fasta = sandbox / (uuid4().hex + '.fasta')
    simulate_fasta(1, 10000, simulated_fasta, include_n=include_n)
    record = next(SeqIO.parse(simulated_fasta, "fasta"))
    sequence = record.seq
    if include_n:
        assert 'N' in sequence
    else:
        assert not 'N' in sequence


def test_index_fasta_with_samtools(nonindexed_fasta):
    samtools_index = local.path(nonindexed_fasta + '.fai')
    assert not samtools_index.exists()
    index_fasta(nonindexed_fasta, 'samtools')
    assert samtools_index.exists()
    samtools_index.delete()


def test_index_fasta_with_bowtie2(nonindexed_fasta):
    bowtie_index = local.path(nonindexed_fasta.with_suffix('.1.bt2'))
    assert not bowtie_index.exists()
    index_fasta(nonindexed_fasta, 'bowtie2')
    assert bowtie_index.exists()
    bowtie_index.delete()


def test_index_fasta_with_all(nonindexed_fasta):
    samtools_index = local.path(nonindexed_fasta + '.fai')
    bowtie_index = local.path(nonindexed_fasta.with_suffix('.1.bt2'))
    assert [not index.exists() for index in [samtools_index, bowtie_index]]
    index_fasta(nonindexed_fasta, 'all')
    assert [index.exists() for index in [samtools_index, bowtie_index]]
