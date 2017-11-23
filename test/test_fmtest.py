"""Test suite for fmbiopy.fmtest"""
def test_cd(cd, startdir):
    assert LocalPath.cwd() != startdir
