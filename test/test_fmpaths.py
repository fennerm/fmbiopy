"""Test fmbiopy.fmpaths"""

import fmbiopy.fmtest as fmtest
import fmbiopy.fmpaths as fmpaths


def test_files_of_type(tmpdir):
    tmpdir = str(tmpdir)
    foo_files = sorted([
        fmtest.gen_tmp(empty=False, suffix='.foo', directory=tmpdir)
        for i in range(0, 4)])
    bar_files = sorted([
        fmtest.gen_tmp(empty=False, suffix='.bar', directory=tmpdir)
        for i in range(0, 4)])

    assert fmpaths.files_of_type(directory=tmpdir, types=['foo']) == foo_files
    assert fmpaths.files_of_type(directory=tmpdir, types=['foo', 'bar']) == \
        sorted(foo_files + bar_files)
