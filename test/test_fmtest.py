"""Test suite for fmbiopy.fmtest"""
from plumbum import local

def test_cd(cd, startdir):
    assert local.cwd != startdir
