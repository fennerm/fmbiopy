"""Test fmbiopy.fmpaths"""

import os

import fmbiopy.fmpaths as fmpaths
from fmbiopy.fmtest import gen_tmp
from fmbiopy.fmtest import nested_dir


def test_contents_to_dict(nested_dir):
    files = ['a.x', 'b.y']
    expected = {}
    for d in fmpaths.listdirs(nested_dir):
        base = os.path.basename(d)
        expected[base] = [os.path.join(d, f) for f in files]

    assert fmpaths.contents_to_dict(nested_dir) == expected


def test_match_files(tmpdir):
    tmpdir = str(tmpdir)
    foo_files = sorted([
        gen_tmp(empty=False, suffix='.foo', directory=tmpdir)
        for i in range(0, 4)])
    bar_files = sorted([
        gen_tmp(empty=False, suffix='.bar', directory=tmpdir)
        for i in range(0, 4)])

    assert fmpaths.match_files(directory=tmpdir, types=['foo']) == foo_files
    assert fmpaths.match_files(directory=tmpdir, types=['foo', 'bar']) == \
        sorted(foo_files + bar_files)
    assert fmpaths.match_files(
            directory=tmpdir, types=['foo', 'bar'],
            substring='.fo') == sorted(foo_files)
