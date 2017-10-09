"""Test fmbiopy.fmpaths"""

import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmtest as fmtest


def test_match_files(tmpdir):
    tmpdir = str(tmpdir)
    foo_files = sorted([
        fmtest.gen_tmp(empty=False, suffix='.foo', directory=tmpdir)
        for i in range(0, 4)])
    bar_files = sorted([
        fmtest.gen_tmp(empty=False, suffix='.bar', directory=tmpdir)
        for i in range(0, 4)])

    assert fmpaths.match_files(directory=tmpdir, types=['foo']) == foo_files
    assert fmpaths.match_files(directory=tmpdir, types=['foo', 'bar']) == \
        sorted(foo_files + bar_files)
    assert fmpaths.match_files(
            directory=tmpdir, types=['foo', 'bar'],
            substring='.fo') == sorted(foo_files)
