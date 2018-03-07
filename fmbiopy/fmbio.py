'''Functions for operating on bioinformatics files'''
from __future__ import print_function

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from numpy.random import choice
from plumbum import (
    FG,
    local,
)
from plumbum.cmd import (
    bowtie2,
    samtools,
    wc,
    zcat,
)


from fmbiopy.fmpaths import is_empty
from fmbiopy.fmsystem import capture_stdout


def align_and_sort(idx, fastq1, fastq2, output_bam, unpaired_fastq=None,
                   max_insert_size=1000, preset='sensitive',
                   threads=1):
    '''Align with bowtie2, convert to bam, sort and index'''
    bowtie2_args = ['-x', idx, '-1', fastq1, '-2', fastq2, '--' + preset,
                    '-p', str(threads), '-X', str(max_insert_size)]
    if unpaired_fastq:
        bowtie2_args += ['-U', unpaired_fastq]

    (bowtie2.__getitem__(bowtie2_args) |
     samtools['view', '-bh', '-'] |
     samtools['sort'] > output_bam) & FG
    samtools['index', output_bam]()


def is_fastq(filename):
    '''Return True if the file is a fastq file'''
    filename = local.path(filename)
    fastq_extensions = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']
    return any([filename.name.lower().endswith(ext) for ext in
                fastq_extensions])


def count_reads(filename):
    """Calculate the number of reads in a fastq, fastq.gz or bam file

    Results written to stdout

    Parameters
    ----------
    filename: Pathlike
    A file path

    Returns
    -------
    int
        The number of reads in the file
    """
    filename = local.path(filename)
    if is_fastq(filename):
        if filename.suffix == '.gz':
            nreads = zcat[filename] | wc['-l']
        else:
            nreads = wc['-l'] < filename
        n = int(int(capture_stdout(nreads)[0]) / 4)
    elif filename.suffix.lower() in ['.sam', '.bam']:
        nreads = samtools['view', '-c', filename]
        n = int(capture_stdout(nreads)[0])
    return n


def index_fasta(filename, method='samtools'):
    '''Index a .fasta file

    Parameters
    ----------
    filename: pathlike
    Path to the fasta file
    method: string, optional
    Indexing method (Possible values: ['samtools', 'bowtie2', 'all'])
    '''
    if method in ['samtools', 'all']:
        samtools['faidx', filename]()
    if method in ['bowtie2', 'all']:
        bowtie_build = local['bowtie2-build']
        bowtie_build[filename, filename.with_suffix('')]()


def merge_bams(bams, output_file, sort_by="index"):
    '''Merge a list of .bam files into a single sorted, indexed .bam file'''

    if sort_by not in ['index', 'name']:
        raise ValueError(
            "sort_by must be one of ['name', 'index'], not %s" % sort_by)

    # Empty bams cause samtools error, exclude them
    bams = [bam for bam in bams if not is_empty(bam)]

    merge_args = ['cat'] + bams
    command = samtools.__getitem__(merge_args)

    if sort_by == 'index':
        command = command | samtools['sort']
    else:
        command = command | samtools['sort', '-n']

    (command > output_file) & FG


def to_fastq(bam, output_prefix, zipped=True):
    '''Convert a bam or sam to .fastq files'''
    reads = {}
    reads['fwd'] = local.path(output_prefix + '.R1.fastq')
    reads['rev'] = local.path(output_prefix + '.R2.fastq')
    reads['unpaired'] = local.path(output_prefix + '.unpaired.fastq')

    bbmap_repair = local['repair.sh']
    bbmap_reformat = local['reformat.sh']

    if zipped:
        for k, v in reads.items():
            reads[k] = local.path(v + '.gz')

    # Piping through bbmap repair and reformat ensures that pairing info is
    # correct
    (samtools['fastq', bam] |
     bbmap_repair['in=stdin.fastq', 'out=stdout', 'outs=' + reads['unpaired']] |
     bbmap_reformat['in=stdin.fastq', 'int=t', 'out1=' + reads['fwd'],
                    'out2=' + reads['rev']]
     ) & FG


def rand_dna_seq(length, base_probs=[0.24, 0.24, 0.24, 0.24, 0.04]):
    '''Generate a random DNA sequence

    Parameters
    ----------
    length: int
        Length of the sequence
    base_probs: List[float]
        Base probabilities in following order [A, C, G, T, N]

    Returns
    -------
    str
    '''
    bases = ['A', 'C', 'G', 'T', 'N']
    sequence = ''.join(choice(bases, length, p=base_probs))
    return sequence


def simulate_fasta(num_sequences, contig_length, output_file, include_n=True):
    '''Simulate a multifasta file

    Parameters
    ----------
    num_sequences: int
        Number of sequences to simulate
    contig_length: int
        Length of each contig
    output_file: pathlike
        Path for the output fasta
    include_n: bool
        If True, 'N' nucleotides are included in the fasta. 'N' nucleotides are
        included rarely(~1 every 100 bases)
    '''
    output_file = local.path(output_file)
    name_root = 'simulated_fasta_' + str(contig_length) + '_'
    if include_n:
        base_probs = [0.2475] * 4 + [0.01]
    else:
        base_probs = [0.25] * 4 + [0]
    with output_file.open('w') as f:
        for i in range(num_sequences):
            sequence = rand_dna_seq(contig_length, base_probs)
            seqid = name_root + str(i)
            record = SeqRecord(Seq(sequence, IUPAC.ambiguous_dna), id=seqid)
            SeqIO.write(record, f, 'fasta')
