"""Test fmbiopy.biofile"""

from pathlib import Path
from typing import (
        Callable,
        Type,
        )

from pytest import (
        fixture,
        raises,
        )

from fmbiopy.biofile import *
from fmbiopy.fmclass import list_classes
from fmbiopy.fmpaths import (
        add_suffix,
        as_paths,
        )


InstanceMap = Callable[[Type[Biofile], str], Biofile]


@fixture(scope='module')
def instance_of(
        example_file: Callable[[str, str], Path]
        )-> InstanceMap:
    def _make_test_instance(
            cls: Type[Biofile],
            size: str)-> Biofile:

        input_example = example_file(cls.input_type, size)

        return cls(input_example)
    return _make_test_instance


@fixture(
        scope='module',
        params=list_classes(
            'biofile',
            package='fmbiopy',
            of_type=['Biofile']))
def biofiles(request):
    return request.param


@fixture(scope='module')
def inst_biofiles(instance_of, biofiles):
    """Return instances of all Biofile types"""
    return instance_of(biofiles, 'tiny')


def test_type_to_class():
    assert type_to_class('fastq') == Fastq
    assert type_to_class('fasta') == Fasta
    assert type_to_class('foo') == Biofile
    assert not type_to_class(None)


class TestBiofile(object):
    """Test `Biofile` class"""

    def test_extension_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'extensions')

    def test_input_type_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'input_type')

    def test_name_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'name')

    def test_path_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'path')

    def test_gzipped_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'gzipped')

    def test_empty_input_raises_value_err(self, biofiles):
        with raises(TypeError):
            biofiles(Path(''))

    def test_undeclared_gzip_raises_gzip_error(self, dat):
        with raises(GzipStatusError):
            Biofile(dat['tiny']['zipped_fwd_reads'][0]).validate()

    def test_list_input_raises_type_error(self, example_file, biofiles):
        with raises(AttributeError):
            biofiles([example_file(biofiles.input_type, 'tiny')])

    def test_if_files_dont_exist_raise_err(self):
        with raises(FileNotFoundError):
            Biofile(Path('i_dont_exist.fa')).validate()

    def test_empty_files_raises_err(self, dat):
        bf = Biofile(dat['tiny']['empty_reads'][0])
        with raises(EmptyFileError):
            bf.validate()

    def test_possibly_empty_prevents_error(self, dat):
        bf = Biofile(dat['tiny']['empty_reads'][0], possibly_empty=True)
        assert bf.validate()

    def test_incorrect_extension_raises_extension_err(self, biofiles, dat):
        read_path = dat['tiny']['fwd_reads'][0]
        incorrect_suffix = add_suffix(read_path, '.foobar')

        if biofiles.extensions != ['ANY']:
            with raises(FileExtensionError):
                biofiles(incorrect_suffix)


@fixture()
def diff_prefix(dat):
    return BiofileGroup(dat['tiny']['diff_prefix'], 'fasta')


@fixture()
def fwd_reads(dat):
    return BiofileGroup(dat['tiny']['fwd_reads'], 'fastq')


@fixture()
def rev_reads(dat):
    return BiofileGroup(dat['tiny']['rev_reads'], 'fastq')


@fixture()
def assemblies(dat):
    return BiofileGroup(dat['tiny']['assemblies'], 'fasta')


class TestBiofileGroup(object):

    def test_empty_input_raises_value_err(self):
        with raises(ValueError):
            BiofileGroup([], filetype='fasta')

    def test_access_returns_path(self, dat, assemblies):
        fasta_paths = dat['tiny']['assemblies']
        expect = fasta_paths[0]
        assert assemblies[0] == expect

    def test_length_method(self, dat, assemblies):
        fasta_paths = dat['tiny']['assemblies']
        actual = len(assemblies)
        expect = len(fasta_paths)
        assert actual == expect

    def test_single_path_input_raises_type_err(self):
        with raises(TypeError):
            BiofileGroup(Path('foo.fa'), filetype='fasta')

    def test_biofilegroups_can_be_zipped(self, assemblies, fwd_reads):
        max_index = max(len(assemblies), len(fwd_reads))
        for fa, reads, i in zip(assemblies, fwd_reads, range(max_index)):
            assert fa == assemblies[i]
            assert reads == fwd_reads[i]

    def test_length_nonexist_doesnt_raise_error(self):
        paths = as_paths(['foo.fa', 'bar.fa'])
        len(BiofileGroup(paths, filetype='fasta'))

    def test_equality_operator(self, assemblies, diff_prefix):
        assert assemblies != diff_prefix

    def test_different_extensions_raises_value_err(self):
        with raises(FileExtensionsNotSameError):
            paths = as_paths(['a.fa', 'b.fasta', 'c.fa'])
            BiofileGroup(paths, filetype='fasta')

    def test_optional_param_are_passed_to_biofile(self, dat):
        with raises(GzipStatusError):
            BiofileGroup(
                    dat['tiny']['fwd_reads'],
                    filetype='fastq',
                    gzipped=True)


class TestMatchedPrefixGroup():
    def test_diff_prefixes_raise_err(self, diff_prefix, fwd_reads, rev_reads):
        with raises(PrefixMatchError):
            MatchedPrefixGroup([diff_prefix, fwd_reads, rev_reads])

    def test_same_filename_raise_value_error(self, fwd_reads):
        with raises(DuplicateFilegroupError):
            MatchedPrefixGroup([fwd_reads, fwd_reads])


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
