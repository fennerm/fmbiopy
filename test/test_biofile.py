from fmbiopy.biofile import BioFileGroup
from fmbiopy.biofile import FastqGroup as Fastq
from fmbiopy.biofile import MatchedPrefixGroup as MatchedPrefixGroup
from fmbiopy.biofile import PairedFastqGroup as PairedFastq
from fmbiopy.fmpaths import add_suffix
import os
import pytest


class TestBioFileGroup(object):

    def test_if_file_exists_then_access_works(self, fasta_paths, fasta):
        assert fasta[1] == fasta_paths[1]

    def test_extension_setter(self, fasta, ):
        assert fasta._extensions == ['fasta'] * 4

    def test_name_setter(self, fasta):
        assert fasta.names == ['FA_SC', 'FB_SC', 'IA_SC', 'IC_SC']

    def test_empty_input_raises_value_err(self):
        with pytest.raises(ValueError):
            BioFileGroup([])

    def test_undeclared_gzip_raises_type_err(self, read_paths):
        with pytest.raises(TypeError):
            BioFileGroup(read_paths[0])

    def test_length_method(self, fasta, fasta_paths):
        assert len(fasta) == len(fasta_paths)

    def test_string_input_raises_type_err(self):
        with pytest.raises(TypeError) as exc:
            BioFileGroup('foo')

    def test_biofiles_can_be_zipped(self, fasta, readfiles):
        max_index = max(len(fasta), len(readfiles))
        for assembly, reads, i in zip(fasta, readfiles, range(max_index)):
            assert assembly == fasta[i]
            assert reads == readfiles[i]

    def test_if_files_dont_exist_raise_os_error(self, nonexistant_fasta):
        with pytest.raises(OSError) as exc:
            nonexistant_fasta[1]

    def test_length_nonexist_doesnt_raise_os_error(self, nonexistant_fasta):
        len(nonexistant_fasta)

    def test_empty_files_raise_value_error(self, empty_paths):
        bfg = BioFileGroup(empty_paths)
        with pytest.raises(ValueError):
            bfg[1]

    def test_possibly_empty_prevents_error(self, empty_paths):
        bfg = BioFileGroup(empty_paths, possibly_empty=True)
        assert bfg[0] == empty_paths[0]

    def test_equality_operator(self, fasta, diff_prefix):
        assert fasta == fasta
        assert fasta != diff_prefix


class TestFastq():
    def test_if_incorrect_extension_raises_value_error(self, read_paths):
        incorrect_suffix = [add_suffix(p, '.x') for p in read_paths[0]]

        fq = Fastq(incorrect_suffix)
        with pytest.raises(ValueError) as exc:

            fq[1]


class TestPairedFastq():

    def test_valid_input_access(self, paired_fastq, fwd_fastq, rev_fastq):
        assert paired_fastq[0] == list([fwd_fastq[0], rev_fastq[0]])

    def test_different_lengths_raises_value_error(self, read_paths):
        fwd_fastq = Fastq(read_paths[0], gzipped=True)
        rev_fastq = Fastq(read_paths[1][1:2], gzipped=True)
        with pytest.raises(ValueError) as exc:
            paired_fastq = PairedFastq(fwd_fastq, rev_fastq)

    def test_length_method(self, fwd_fastq, paired_fastq):
        assert len(paired_fastq) == len(fwd_fastq)


class TestBowtie2Index():
    def test_bowtie2_index_name_setter(self, bowtie2_indices):
        testdir = os.path.abspath('sandbox/bowtie2_indices/')
        assert bowtie2_indices._index_prefixes == [
                testdir+'/FA_SC.scaffolds',
                testdir+'/FB_SC.scaffolds',
                testdir+'/IA_SC.scaffolds',
                testdir+'/IC_SC.scaffolds'
                ]


class TestIndexedFasta():
    def test_attributes(self, fasta, samtools_indices, bowtie2_indices,
            indexed_fasta):
        assert indexed_fasta.fasta == fasta
        assert indexed_fasta.indices == tuple([samtools_indices,
                bowtie2_indices])

    def test_access(self, fasta, samtools_indices, bowtie2_indices,
            indexed_fasta):
        assert indexed_fasta[0] == list([fasta[0], samtools_indices[0],
                bowtie2_indices[0]])


class TestMatchedPrefixGroup():
    def test_diff_prefixes_raise_value_error(self, diff_prefix, fwd_fastq,
            rev_fastq):
        with pytest.raises(ValueError):
            MatchedPrefixGroup([diff_prefix, fwd_fastq, rev_fastq])

    def test_same_filename_raise_value_error(self, fwd_fastq):
        with pytest.raises(ValueError):
            MatchedPrefixGroup([fwd_fastq, fwd_fastq])
