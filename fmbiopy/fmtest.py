"""Set of functions to aid in testing

Modules which import fmtest must also import load_sandbox explicitely.
"""

import glob
import os
import pytest
import shutil
import tempfile
from typing import Dict
from typing import Generator
from typing import Iterator
from typing import List
from typing import Tuple

import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmsystem as fmsystem


def gen_tmp(
        empty: bool = True,
        suffix: str = '',
        directory: str = 'sandbox') -> str:
    """Generate a named temporary file.

    Warning: These files need to be deleted manually if a non-temporary
    directory is used.

    Parameters
    ----------
    empty
        If True, the file is empty. Otherwise it has content.
    suffix, optional
        If defined, the generated files will have the given extension
    directory, optional
        If defined, the generated files will be produced in the given
        directory. By default the directory produced by load_sandbox is used.

    Returns
    -------
    The path to the created temporary file
    """

    tmpfile = tempfile.NamedTemporaryFile(
            delete=False, dir=directory, suffix=suffix).name

    if not empty:
        with open(tmpfile, 'w') as f:
            f.write('foo')
    else:
        with open(tmpfile, 'w') as f:
            f.write('')
    return tmpfile


def gen_mixed_tmpfiles(*args, **kwargs) -> List[str]:
    """Generate a list of two tempfiles - the first is nonempty

    All arguments are passed to `gen_tmp`
    """
    tmps = [gen_tmp(empty=False, *args, **kwargs)] + \
        [gen_tmp(empty=True, *args, **kwargs)]
    with open(tmps[0], "w") as f:
        f.write("foo")
    return tmps


@pytest.fixture(scope='session', autouse=True)
def load_sandbox() -> Generator:
    """Copy all test data files to the sandbox for the testing session"""
    if os.path.exists('sandbox'):
        shutil.rmtree('sandbox')

    shutil.copytree('testdat', 'sandbox')
    yield
    shutil.rmtree('sandbox')


@pytest.fixture(scope='class', autouse=True)
def initial_test_state() -> Iterator[Tuple[str, List[str], List[str]]]:
    """Stores the initial state of the test data directory"""
    return os.walk('sandbox')


def get_dat() -> Dict[str, List[str]]:
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
