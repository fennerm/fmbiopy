"""Set of functions to aid in testing"""

import glob
import os
import pytest
import shutil
import tempfile
import typing

import fmbiopy.fmpaths as fmpaths


def gen_tmp(empty=True)-> str:
    """Generate a named temporary file

    Parameters
    ----------
    empty
        If True, the file is empty. Otherwise it has content.

    Returns
    -------
    The path to the created temporary file
    """
    tmpfile = tempfile.NamedTemporaryFile(delete=False).name
    if not empty:
        with open(tmpfile, 'w') as f:
            f.write('foo')

    return tmpfile


def gen_mixed_tmpfiles():
    """Generate a list of two tempfiles - the first is nonempty"""
    tmps = [tempfile.NamedTemporaryFile(delete=False).name
            for i in range(0, 2)]
    with open(tmps[0], "w") as f:
        f.write("foo")
    return tmps


@pytest.fixture(scope='session', autouse=True)
def load_sandbox() -> None:
    """Copy all test data files to the sandbox for the testing session"""
    if os.path.exists('sandbox'):
        shutil.rmtree('sandbox')

    shutil.copytree('testdat', 'sandbox')
    yield
    shutil.rmtree('sandbox')


@pytest.fixture(scope='class', autouse=True)
def initial_test_state() -> typing.List[str]:
    """Stores the initial state of the test data directory"""
    return os.walk('sandbox')


def get_dat() -> typing.Dict[str, typing.List[str]]:
    """Create a dictionary of test data

    Assumes test directory is structured such that all test data is stored in
    the test/testdat/sandbox directory. test/testdat can contain any number of
    directories which each store a certain group of data files. `get_dat`
    represents this structure as a dictionary with subdirectories of sandbox as
    keys and datafile paths as values

    Returns
    -------
    A dictionary of the form Dict[Subdirectories of testdat, files in
    subdirectory].

    Designed to be run from the test directory, which contains a testdat
    directory.
    """

    testdirs = [
            os.path.abspath(d) for d in fmpaths.listdirs('sandbox')]
    dat = {}
    for d in testdirs:
        base = os.path.basename(d)
        dat[base] = sorted(glob.glob(d + '/*'))
    return dat
