import pytest
from fmbiopy.fmpaths import abs_paths, add_suffix
from get_dat import get_dat
from fmbiopy.biofile import (
        BioFileGroup,
        Fastq,
        PairedFastq,
        Bam,
        Sam,
        )

@pytest.fixture
def assembly_paths():
    dat = get_dat()['assemblies']
    return dat

@pytest.fixture
def read_paths():
    dat = get_dat()
    return (dat['fwd_reads'], dat['rev_reads'])

@pytest.fixture
def empty_paths():
    dat = get_dat()
    return (dat['empty'])


class TestBioFileGroup():
    @pytest.fixture
    def assfiles(self, assembly_paths):
        return BioFileGroup(assembly_paths)

    @pytest.fixture
    def readfiles(self, read_paths):
        return BioFileGroup(read_paths[0], gzipped=True)

    def test_properties_set(self, assfiles):
        assfiles.paths
        assfiles.gzipped
        assfiles._validated
        assfiles._accepted_extensions
        assfiles.extensions

    def test_if_file_exists_then_access_works(self, assembly_paths, assfiles):
        assert assfiles.paths[1] == assembly_paths[1]

    def test_empty_input_raises_value_err(self):
        with pytest.raises(ValueError) as exc:
            BioFileGroup([])

    def test_undeclared_gzip_raises_type_err(self, read_paths):
        with pytest.raises(TypeError) as exc:
            BioFileGroup(read_paths[0])

    def test_length_method(self, assfiles, assembly_paths):
        assert len(assfiles) == len(assembly_paths)

    def test_string_input_raises_type_err(self):
        with pytest.raises(TypeError) as exc:
            BioFileGroup('foo')

    def test_biofiles_can_be_zipped(self, assfiles, readfiles):
        max_index = max(len(assfiles), len(readfiles))
        for assembly, reads, i in zip(assfiles, readfiles, range(max_index)):
            assert assembly == assfiles[i]
            assert reads == readfiles[i]

    def test_if_files_dont_exist_raise_os_error(self, assembly_paths):
        nonexistant = ['foo/' + path for path in assembly_paths[0]]

        bfg = BioFileGroup(nonexistant)
        with pytest.raises(OSError) as exc:
            bfg[1]

    def test_empty_files_raise_value_error(self, empty_paths):
        bfg = BioFileGroup(empty_paths)
        with pytest.raises(ValueError) as exc:
            bfg[1]

    def test_possibly_empty_prevents_error(self, empty_paths):
        bfg = BioFileGroup(empty_paths, possibly_empty=True)
        assert bfg[0] == empty_paths[0]


class TestFastq():

    @pytest.fixture
    def reads(self, read_paths):
        return Fastq(read_paths[0], gzipped=True)

    def test_set_attributes(self, reads):
        assert reads._accepted_extensions == ['fastq', 'fq']

    def test_if_incorrect_extension_raises_value_error(self, read_paths):
        incorrect_suffix = add_suffix(read_paths[0], '.x')

        fq = Fastq(incorrect_suffix)
        with pytest.raises(ValueError) as exc:
            fq[1]

