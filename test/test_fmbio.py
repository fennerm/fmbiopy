from __future__ import division
from collections import Counter
from uuid import uuid4

from Bio import SeqIO
from pytest import (
    approx,
    mark,
    warns,
)

from fmbiopy.fmbio import *
from fmbiopy.fmlist import none
from fmbiopy.fmpaths import (
    all_exist,
    apply_is_empty,
    is_empty,
)


@mark.parametrize("with_unpaired", [True, False])
def test_align_and_sort_no_unpaired(sandbox, partial_fasta,
                                    paired_trimmed_fastq, with_unpaired):
    output_bam = sandbox / (uuid4().hex + '.bam')
    bam_index = local.path(output_bam + '.bai')
    args = {
        'idx': partial_fasta['bt2'],
        'fastq1': paired_trimmed_fastq['fwd'],
        'fastq2': paired_trimmed_fastq['rev'],
        'preset': 'very-fast',
        'output_bam': output_bam
    }

    if with_unpaired:
        args['unpaired_fastq'] = paired_trimmed_fastq['unpaired']

    align_and_sort(**args)
    assert not is_empty(output_bam)
    assert not is_empty(bam_index)


def test_align_and_sort_with_unpaired(sandbox, paired_trimmed_fastq):
    pass


def test_count_reads_output_correct_for_bam(untrimmed_bam):
    n = count_reads(untrimmed_bam)
    assert n > 830 and n < 860


def test_count_reads_output_correct_for_empty_bam(empty_bam):
    assert count_reads(empty_bam) == 0


def test_count_reads_output_correct_for_sam(untrimmed_sam):
    n = count_reads(untrimmed_sam)
    assert n > 830 and n < 850


def test_count_reads_output_correct_for_fastq(simulated_reads):
    n = count_reads(simulated_reads['fwd'])
    assert n > 415 and n < 425


@mark.parametrize('sort_by', ['index', 'name'])
def test_merge_bams(mini_bams, sandbox, sort_by):
    output_file = sandbox / (uuid4().hex + '.bam')
    merge_bams(mini_bams, output_file, sort_by=sort_by)
    num_input_reads = sum([count_reads(bam) for bam in mini_bams])
    assert count_reads(output_file) == num_input_reads


def test_merge_bams_with_empty_bam(sandbox, mini_bams, empty_bam):
    output_file = sandbox / (uuid4().hex + '.bam')
    merge_bams(mini_bams + [empty_bam], output_file)


def test_to_fastq(bam_with_orphans, sandbox):
    output_prefix = sandbox / uuid4().hex
    suffixes = ['.R1.fastq.gz', '.R2.fastq.gz', '.unpaired.fastq.gz']
    expected_output = [local.path(output_prefix + suff) for suff in suffixes]
    to_fastq(bam_with_orphans, output_prefix)
    assert all_exist(expected_output) and none(apply_is_empty(expected_output))
    assert count_reads(expected_output[0]) == count_reads(expected_output[1])


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


@mark.parametrize("length", [0, 1, 10])
def test_rand_dna_seq_length(length):
    assert len(rand_dna_seq(length)) == length


def test_rand_dna_seq_contains_expected_nucleotides():
    assert set('ACGTN') == set(rand_dna_seq(10000))


def test_rand_dna_seq_base_probs():
    length = 10000
    base_prob = 0.2
    base_probs = [base_prob] * 5
    seq = rand_dna_seq(length, base_probs)
    base_counts = Counter(seq)
    base_freqs = [v / length for v in base_counts.values()]

    for freq in base_freqs:
        assert freq == approx(base_prob, 0.05)
