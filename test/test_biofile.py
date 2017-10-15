"""Test fmbiopy.biofile"""

from pathlib import Path
import pytest
from typing import Callable
from typing import Type

import fmbiopy.biofile as biofile
import fmbiopy.fmclass as fmclass
import fmbiopy.fmpaths as fmpaths


InstanceMap = Callable[[Type[biofile.Biofile], str], biofile.Biofile]


@pytest.fixture
def instance_of(
        example_file: Callable[[str, str], Path]
        )-> InstanceMap:
    def _make_test_instance(
            cls: Type[biofile.Biofile],
            size: str)-> biofile.Biofile:

        input_example = example_file(cls.input_type, size)

        return cls(input_example)
    return _make_test_instance


@pytest.fixture(
        params=fmclass.list_classes(
            'fmbiopy.biofile',
            of_type='Biofile'))
def biofiles(request):
    return request.param


@pytest.fixture
def inst_biofiles(instance_of, biofiles):
    """Return instances of all Biofile types"""
    return instance_of(biofiles, 'tiny')


# @pytest.fixture(
#         params=fmclass.list_classes(
#             'fmbiopy.biofile',
#             of_type='Biofile'))
# def groups(request, instance_of):
#     biofile_class = request.param
#     return instance_of(biofile_class)


def test_type_to_class():
    assert biofile.type_to_class('fastq') == biofile.Fastq
    assert biofile.type_to_class('fasta') == biofile.Fasta
    assert biofile.type_to_class('foo') == biofile.Biofile
    assert not biofile.type_to_class(None)


class TestBiofile(object):
    """Test `Biofile` class"""

    def test_extension_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'accepted_extensions')

    def test_input_type_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'input_type')

    def test_name_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'name')

    def test_path_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'path')

    def test_gzipped_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'gzipped')

    def test_empty_input_raises_value_err(self, biofiles):
        with pytest.raises(TypeError):
            biofiles(Path(''))

    def test_undeclared_gzip_raises_gzip_error(self, dat):
        with pytest.raises(biofile.GzipStatusError):
            biofile.Biofile(dat['tiny']['zipped_fwd_reads'][0]).validate()

    def test_list_input_raises_type_error(self, example_file, biofiles):
        with pytest.raises(AttributeError):
            biofiles([example_file(biofiles.input_type, 'tiny')])

    def test_if_files_dont_exist_raise_err(self):
        with pytest.raises(FileNotFoundError):
            biofile.Biofile(Path('i_dont_exist.fa')).validate()

    def test_empty_files_raises_err(self, dat):
        bf = biofile.Biofile(dat['tiny']['empty_reads'][0])
        with pytest.raises(biofile.EmptyFileError):
            bf.validate()

    def test_possibly_empty_prevents_error(self, dat):
        bf = biofile.Biofile(dat['tiny']['empty_reads'][0], possibly_empty=True)
        assert bf.validate()

    def test_incorrect_extension_raises_extension_err(self, biofiles, dat):
        read_path = dat['tiny']['fwd_reads'][0]
        incorrect_suffix = fmpaths.add_suffix(read_path, '.foobar')

        if biofiles.accepted_extensions != ['ANY']:
            with pytest.raises(biofile.FileExtensionError):
                biofiles(incorrect_suffix)

# class TestBiofileGroup(object):
#
#     def test_empty_input_raises_value_err(self):
#         with pytest.raises(ValueError):
#             biofile.BiofileGroup([])
#
#     def test_access_returns_path(self, dat):
#         fasta_paths = dat['assemblies']
#         assert biofile.BiofileGroup(fasta_paths)[0] == fasta_paths[0]
#
#     def test_length_method(self, dat):
#         fasta_paths = dat['assemblies']
#         assert len(biofile.Biofilegroup(fasta_paths)) == len(fasta_paths)
#
#     def test_string_input_raises_type_err(self):
#         with pytest.raises(TypeError):
#             biofile.BiofileGroup('foo')
#
#     def test_biofilegroups_can_be_zipped(self, dat):
#         fasta_group = biofile.BiofileGroup(dat['assemblies'])
#         fwd_read_group = biofile.BiofileGroup(dat['fwd_reads'])
#         max_index = max(len(fasta_group), len(fwd_read_group))
#         for fa, reads, i in zip(fasta_group, fwd_read_group, range(max_index)):
#             assert fa == fasta_group[i]
#             assert reads == fwd_read_group[i]
#
#     def test_length_nonexist_doesnt_raise_error(self):
#         len(biofile.BiofileGroup(['foo.foo', 'bar.bar']))
#
#     def test_equality_operator(self, dat):
#         a = biofile.BiofileGroup(dat['assemblies'])
#         b = biofile.BiofileGroup(dat['assemblies'])
#         c = biofile.BiofileGroup(dat['diff_prefix'])
#         assert a == b
#         assert a != c
#
#     def test_different_extensions_raises_value_err(self):
#         with pytest.raises(ValueError):
#             biofile.BiofileGroup(['a.fa', 'b.fasta', 'c.fa'])
#
#     def test_one_nonexistant_file_caught_by_validation(self):
#         assert False
#
#     def test_biofile_initialization(self):
#         assert False
#
# class TestMatchedPrefixGroup():
#     def test_diff_prefixes_raise_err(self, diff_prefix, fwd_fastq,
#                                      rev_fastq):
#         with pytest.raises(biofile.PrefixMatchError):
#             biofile.MatchedPrefixGroup([diff_prefix, fwd_fastq, rev_fastq])
#
#     def test_same_filename_raise_value_error(self, fwd_fastq):
#         with pytest.raises(biofile.DuplicateFilegroupError):
#             biofile.MatchedPrefixGroup([fwd_fastq, fwd_fastq])
#
#     def test_length_method(self, fwd_fastq, paired_fastq):
#         assert len(paired_fastq) == len(fwd_fastq)
#
#
# class TestBowtie2Index():
#     def test_bowtie2_index_name_setter(self, bowtie2_indices):
#         testdir = os.path.abspath('sandbox/bowtie2_indices/')
#         assert bowtie2_indices._index_prefixes == [
#                 testdir + '/FA_SC.scaffolds',
#                 testdir + '/FB_SC.scaffolds',
#                 testdir + '/IA_SC.scaffolds',
#                 testdir + '/IC_SC.scaffolds']
#
#
# class TestIndexedFasta():
#     def test_attributes(self, fasta, samtools_indices, bowtie2_indices,
#                         indexed_fasta):
#         assert indexed_fasta.fasta == fasta
#         assert indexed_fasta.indices == tuple([samtools_indices,
#                                                bowtie2_indices])
#
#     def test_access(
#             self, fasta, samtools_indices, bowtie2_indices, indexed_fasta):
#         assert indexed_fasta[0] == [fasta[0], samtools_indices[0],
#                                     bowtie2_indices[0]]
