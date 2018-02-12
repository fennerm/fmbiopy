'''Functions for operating on bioinformatics files'''
from __future__ import print_function
from numpy.random import choice
from uuid import uuid4

from plumbum import (
    FG,
    local,
)
from plumbum.cmd import (
    samtools,
    wc,
    zcat,
)

from fmbiopy.fmpaths import delete
from fmbiopy.fmsystem import capture_stdout


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
    if len(bams) > 1:
        merge_args = ['cat'] + bams

        tmpfile = local.path(uuid4().hex + '.bam')
        with delete(tmpfile):
            (samtools.__getitem__(merge_args) | samtools['sort', '-n'] |
             samtools['fixmate', '-', tmpfile]) & FG

            if sort_by == 'index':
                (samtools['sort', tmpfile] > output_file) & FG
            elif sort_by == 'name':
                tmpfile.move(output_file)
    elif len(bams) == 1:
        (samtools['sort', '-n', bams[0]] |
         samtools['fixmate', '-', output_file]) & FG
    else:
        raise ValueError("len(bams) must be > 0")


def to_fastq(bam, output_prefix, zipped=True):
    '''Convert a bam or sam to .fastq files'''
    reads = {}
    reads['fwd'] = local.path(output_prefix + '.R1.fastq')
    reads['rev'] = local.path(output_prefix + '.R2.fastq')
    reads['unpaired'] = local.path(output_prefix + '.unpaired.fastq')

    if zipped:
        for k, v in reads.items():
            reads[k] = local.path(v + '.gz')
    samtools['fastq', '-O', '-s', reads['unpaired'], '-1', reads['fwd'], '-2',
             reads['rev'], bam] & FG


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
    included rarely (~1 every 100 bases)
    '''
    output_file = local.path(output_file)
    name_base = '>simulated_fasta_' + str(contig_length) + '_'
    if include_n:
        bases = ['A', 'C', 'G', 'T', 'N']
        base_probs = [0.2475] * 4 + [0.01]
    else:
        bases = ['A', 'C', 'G', 'T']
        base_probs = [0.25] * 4
    with output_file.open('w') as f:
        for i in range(num_sequences):
            sequence = ''.join(choice(bases, contig_length, p=base_probs))
            print(name_base + str(i), file=f)
            print(sequence, file=f)
