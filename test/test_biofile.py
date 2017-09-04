import pytest
from fmbiopy.fmpaths import abs_paths
from get_dat import get_dat
from fmbiopy.biofile import (
        BioFileGroup,
        FastqFileGroup,
        BamFileGroup,
        SamFileGroup,
        )

@pytest.fixture
def assembly_paths():
    dat = get_dat()['assemblies']
    return dat

@pytest.fixture
def read_paths():
    dat = get_dat()
    return (dat['fwd_reads'], dat['rev_reads'])


class TestBioFileGroup():
    @pytest.fixture
    def biofiles_exist(self, assembly_paths):
        return BioFileGroup(assembly_paths)

    def test_properties_set(self, biofiles_exist):
        biofiles_exist.paths
        biofiles_exist.gzipped
        biofiles_exist._validated

    def test_if_file_exists_then_access_works(self, assembly_paths,
            biofiles_exist):
        assert biofiles_exist.paths[1] == assembly_paths[1]

    def test_empty_input_raises_value_err(self):
        with pytest.raises(ValueError) as exc:
            BioFileGroup([])

    def test_undeclared_gzip_raises_type_err(self, read_paths):
        with pytest.raises(TypeError) as exc:
            BioFileGroup(read_paths[0])

    def test_length_method(self, biofiles_exist, assembly_paths):
        assert biofiles_exist.length() == len(assembly_paths)

    def test_string_input_raises_type_err(self):
        with pytest.raises(TypeError) as exc:
            BioFileGroup('foo')
