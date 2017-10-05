import pytest
from typing import List
from typing import Tuple

from fmbiopy.biofile import BioFileGroup
from fmbiopy.biofile import Bowtie2IndexGroup as Bowtie2Index
from fmbiopy.biofile import FastaGroup as Fasta
from fmbiopy.biofile import FastqGroup as Fastq
from fmbiopy.biofile import IndexedFastaGroup as IndexedFasta
from fmbiopy.biofile import PairedFastqGroup as PairedFastq
from fmbiopy.biofile import SamtoolsFAIndexGroup as SamtoolsFAIndex
import fmbiopy.fmtest as fmtest
from fmbiopy.fmtest import get_dat
from fmbiopy.fmtest import initial_test_state
from fmbiopy.fmtest import load_sandbox


@pytest.fixture
def fasta_paths() -> List[str]:
    dat = fmtest.get_dat()['assemblies']
    return dat


@pytest.fixture
def read_paths() -> Tuple[List[str], List[str]]:
    dat = fmtest.get_dat()
    return (dat['fwd_reads'], dat['rev_reads'])


@pytest.fixture
def diff_prefix_paths():
    return fmtest.get_dat()['diff_prefix']


@pytest.fixture
def diff_prefix(diff_prefix_paths):
    return Fasta(diff_prefix_paths)


@pytest.fixture
def empty_paths():
    return fmtest.get_dat()['empty']


@pytest.fixture
def fasta(fasta_paths):
    return BioFileGroup(fasta_paths)


@pytest.fixture
def bowtie_index_paths():
    return fmtest.get_dat()['bowtie2_indices']


@pytest.fixture
def samtools_index_paths():
    return fmtest.get_dat()['faindices']


@pytest.fixture
def fwd_fastq(read_paths):
    return Fastq(read_paths[0], gzipped=True)


@pytest.fixture
def rev_fastq(read_paths):
    return Fastq(read_paths[1], gzipped=True)


@pytest.fixture
def bowtie2_indices(bowtie_index_paths):
    return Bowtie2Index(bowtie_index_paths)


@pytest.fixture
def samtools_indices(samtools_index_paths):
    return SamtoolsFAIndex(samtools_index_paths)


@pytest.fixture
def paired_fastq(read_paths, fwd_fastq, rev_fastq):
    return PairedFastq(fwd_fastq, rev_fastq)


@pytest.fixture
def nonexistant_fasta(fasta_paths):
    nonexistant = ['foo/' + path for path in fasta_paths[0]]

    return BioFileGroup(nonexistant)


@pytest.fixture
def readfiles(read_paths):
    return BioFileGroup(read_paths[0], gzipped=True)


@pytest.fixture
def indexed_fasta(fasta, samtools_indices, bowtie2_indices):
    return IndexedFasta(fasta, samtools_indices, bowtie2_indices)
