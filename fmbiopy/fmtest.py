"""Set of functions to aid in testing

Modules which import must also import load_sandbox explicitely.
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

import fmbiopy.fmlist as fmlist
import fmbiopy.fmpaths as fmpaths


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

    def _ignore_git(*args):
        """Copying the git directory is unnecessary so we ignore it

        This function is used by `shutil.copytree`"""
        return '.git'

    if os.path.exists('sandbox'):
        shutil.rmtree('sandbox')

    shutil.copytree('testdat', 'sandbox', ignore=_ignore_git)
    yield
    if os.path.exists('sandbox'):
        shutil.rmtree('sandbox')


@pytest.fixture(scope='class', autouse=True)
def initial_test_state() -> Iterator[Tuple[str, List[str], List[str]]]:
    """Stores the initial state of the test data directory"""
    return os.walk('sandbox')


@pytest.fixture
def example_file(dat, tmpdir):
    """Given a file extension, return a testfile of that type"""

    def get_example_file(filetype):
        if filetype == 'fasta':
            return dat['assemblies'][0]
        elif filetype == ('fastq', 'fastq'):
            return (dat['fwd_reads'][0], dat['rev_reads'][0])
        elif filetype == 'fastq':
            return dat['fwd_reads'][0]
        elif filetype == 'fai':
            return dat['faindices'][0]
        # elif filetype == 'sam':
        #     return dat['sam'][0]
        # elif filetype == 'bam':
        #     return dat['bam'][0]
        elif filetype == 'gz':
            return dat['zipped_fwd_reads'][0]
        return gen_tmp(empty=False, directory=tmpdir, suffix='.foo')

    return get_example_file


@pytest.fixture
def instance_of(example_file):
    """Given a class name, return an instance of the task

    Only works for classes with the class attributes input_type and/or
    output_type. Extra parameters cannot be passed. Designed for testing groups
    of closely related classes which are all initialized using the same basic
    process but with different types of input files.
    """
    def make_test_instance(class_name):
        input_example = [example_file(t) for t in class_name.input_type]
        input_example = fmlist.flatten(input_example)

        try:
            if len(input_example) == 1 and class_name.output_type == ['']:
                output_example = [fmpaths.remove_suffix(input_example[0])]
            else:
                # If we have multiple inputs, the output suffix is added to
                # the first input as in ruffus
                input_prefix = fmpaths.remove_suffix(input_example[0])
                output_example = []
                for typ in class_name.output_type:
                    output_example.append(
                            fmpaths.add_suffix(input_prefix, '.' + typ))

            if len(input_example) == 1:
                input_example = input_example[0]
            if len(output_example) == 1:
                output_example = output_example[0]
            return class_name(input_example, output_example)
        except AttributeError:
            if len(input_example) == 1:
                input_example = input_example[0]
            return class_name(input_example)
    return make_test_instance


@pytest.fixture(scope='session')
def dat() -> Dict[str, List[str]]:
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
