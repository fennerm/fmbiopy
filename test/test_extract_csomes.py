from __future__ import print_function
from uuid import uuid4

from Bio import SeqIO
from plumbum.cmd import (
    bowtie2,
    samtools,
)
from pytest import (
    fixture,
    mark,
)

from fmbiopy.fmbio import count_reads
from fmbiopy.fmtest import assert_script_produces_files


@fixture
def include_list(sandbox, partial_fasta):
    output_file = sandbox / (uuid4().hex + '.txt')
    with output_file.open('w') as f:
        contigs = SeqIO.parse(partial_fasta['fasta'], 'fasta')
        for i, contig in enumerate(contigs):
            if i > 2 and i < 6:
                print(contig.id, file=f)
    return output_file


@mark.parametrize('output_format', ['fastq', 'bam'])
def test_extract_csomes(sandbox, partial_fasta, trimmed_bam, output_format,
                        include_list):
    output_prefix = sandbox / uuid4().hex

    if output_format == 'fastq':
        suffixes = ['.R1.fastq.gz', '.R2.fastq.gz', '.unpaired.fastq.gz']
    elif output_format == 'bam':
        suffixes = ['.bam']

    expected_output = [output_prefix + suff for suff in suffixes]

    assert_script_produces_files("bin/extract_csomes.py",
                                 args=["-o", output_prefix, '-f', output_format,
                                       '-i', include_list, '-p', '2',
                                       trimmed_bam],
                                 output=expected_output)

    num_input_reads = count_reads(trimmed_bam)
    num_output_reads = sum([count_reads(f) for f in expected_output])
    num_output_reads = [count_reads(f) for f in expected_output]

    if output_format == 'fastq':
        assert num_output_reads[0] == num_output_reads[1]
        # Check that files are mostly paired
        assert num_output_reads[0] > num_output_reads[2]

    assert sum(num_output_reads) < num_input_reads
