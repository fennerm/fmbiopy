"""Test suite for fmbiopy.fmtest"""
from pathlib import Path

def test_cd(cd, startdir):
    assert Path.cwd() != startdir
