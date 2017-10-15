"""Test fmbiopy.fmpaths"""

from pathlib import Path
import os
import pytest

import fmbiopy.fmpaths as fmpaths
from fmbiopy.fmtest import full_dir
from fmbiopy.fmtest import gen_tmp
from fmbiopy.fmtest import nested_dir
from fmbiopy.fmtest import tmpdir


def test_prefix():
    assert fmpaths.prefix(Path('foo/bar/name.a.b')) == 'name'


def test_as_dict():
    # Cant figure out a way to test this one without just rewriting the
    # function
    pass


@pytest.fixture()
def bowtie2_suffixes():
    return list(['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
                 '.rev.2.bt2'])


def test_get_bowtie2_indices(bowtie2_suffixes):
    prefix = 'foo'
    expected_paths = [prefix + suffix for suffix in bowtie2_suffixes]
    dot_paths = fmpaths.get_bowtie2_indices(prefix)
    actual_paths = [str(path) for path in dot_paths]
    assert actual_paths == expected_paths

@pytest.mark.parametrize("ext,substr,expect", [
        (None, None, [Path('a.x'), Path('b.y'), Path('b.x')]),
        ('y', None, [Path('b.y')]),
        (None, 'a', [Path('a.x')]),
        (None, 'k', []),
        ('x', 'b', [Path('b.x')])])
def test_find(ext, substr, expect, full_dir):
    if expect:
        expect = sorted([full_dir.joinpath(exp) for exp in expect])
    assert fmpaths.find(full_dir, extensions=ext, substring=substr) == expect
