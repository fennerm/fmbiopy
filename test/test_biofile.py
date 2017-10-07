import os
import pytest

import fmbiopy.biofile as biofile
import fmbiopy.fmpaths as fmpaths


class TestBioFileGroup(object):

    def test_if_file_exists_then_access_works(self, fasta_paths, fasta):
        assert fasta[1] == fasta_paths[1]

    def test_extension_setter(self, fasta, ):
        assert fasta._extensions == ['fasta'] * 4

    def test_name_setter(self, fasta):
        assert fasta.names == ['FA_SC', 'FB_SC', 'IA_SC', 'IC_SC']

    def test_empty_input_raises_value_err(self):
        with pytest.raises(ValueError):
            biofile.BioFileGroup([])

    def test_undeclared_gzip_raises_gzip_error(self, read_paths):
        with pytest.raises(biofile.GzipStatusError):
            biofile.BioFileGroup(read_paths[0])

    def test_length_method(self, fasta, fasta_paths):
        assert len(fasta) == len(fasta_paths)

    def test_string_input_raises_type_err(self):
        with pytest.raises(TypeError):
            biofile.BioFileGroup('foo')

    def test_biofiles_can_be_zipped(self, fasta, readfiles):
        max_index = max(len(fasta), len(readfiles))
        for assembly, reads, i in zip(fasta, readfiles, range(max_index)):
            assert assembly == fasta[i]
            assert reads == readfiles[i]

    def test_if_files_dont_exist_raise_os_error(self, nonexistant_fasta):
        with pytest.raises(OSError):
            nonexistant_fasta[1]

    def test_length_nonexist_doesnt_raise_os_error(self, nonexistant_fasta):
        len(nonexistant_fasta)

    def test_empty_files_raises_err(self, empty_paths):
        bfg = biofile.BioFileGroup(empty_paths)
        with pytest.raises(biofile.EmptyFileError):
            bfg[1]

    def test_possibly_empty_prevents_error(self, empty_paths):
        bfg = biofile.BioFileGroup(empty_paths, possibly_empty=True)
        assert bfg[0] == empty_paths[0]

    def test_equality_operator(self, fasta, diff_prefix):
        assert fasta == fasta
        assert fasta != diff_prefix


class TestFastq():
    def test_if_incorrect_extension_raises_extension_err(self, read_paths):
        incorrect_suffix = [fmpaths.add_suffix(p, '.x') for p in read_paths[0]]

        fq = biofile.FastqGroup(incorrect_suffix)
        with pytest.raises(biofile.FileExtensionError):
            fq[1]


class TestPairedFastq():

    def test_valid_input_access(self, paired_fastq, fwd_fastq, rev_fastq):
        assert paired_fastq[0] == list([fwd_fastq[0], rev_fastq[0]])

    def test_different_lengths_raises_err(self, read_paths):
        fwd_fastq = biofile.FastqGroup([read_paths[0][0]] * 3, gzipped=True)
        rev_fastq = biofile.FastqGroup([read_paths[1][0]] * 2, gzipped=True)
        with pytest.raises(biofile.GroupLengthError):
            biofile.PairedFastqGroup(fwd_fastq, rev_fastq)

    def test_length_method(self, fwd_fastq, paired_fastq):
        assert len(paired_fastq) == len(fwd_fastq)


class TestBowtie2Index():
    def test_bowtie2_index_name_setter(self, bowtie2_indices):
        testdir = os.path.abspath('sandbox/bowtie2_indices/')
        assert bowtie2_indices._index_prefixes == [
                testdir + '/FA_SC.scaffolds',
                testdir + '/FB_SC.scaffolds',
                testdir + '/IA_SC.scaffolds',
                testdir + '/IC_SC.scaffolds']


class TestIndexedFasta():
    def test_attributes(self, fasta, samtools_indices, bowtie2_indices,
                        indexed_fasta):
        assert indexed_fasta.fasta == fasta
        assert indexed_fasta.indices == tuple([samtools_indices,
                                               bowtie2_indices])

    def test_access(
            self, fasta, samtools_indices, bowtie2_indices, indexed_fasta):
        assert indexed_fasta[0] == [fasta[0], samtools_indices[0],
                                    bowtie2_indices[0]]


class TestMatchedPrefixGroup():
    def test_diff_prefixes_raise_err(self, diff_prefix, fwd_fastq,
                                     rev_fastq):
        with pytest.raises(biofile.PrefixMatchError):
            biofile.MatchedPrefixGroup([diff_prefix, fwd_fastq, rev_fastq])

    def test_same_filename_raise_value_error(self, fwd_fastq):
        with pytest.raises(biofile.DuplicateFilegroupError):
            biofile.MatchedPrefixGroup([fwd_fastq, fwd_fastq])
