"""Set of functions to aid in testing

Pytest functions must be imported explicitely.
"""

import os
import shutil
import tempfile
from typing import Callable
from typing import Dict
from typing import Generator
from typing import Iterator
from typing import List
from typing import Tuple
from typing import Type
from typing import TypeVar

import py
import pytest

import fmbiopy.biofile as biofile
import fmbiopy.fmlist as fmlist
import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmruffus as fmruffus


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


@pytest.fixture(scope='session', autouse=True)
def load_sandbox() -> Generator:
    """Copy all test data files to the sandbox for the testing session

    Yields
    ------
    A temporary directory with the contents copied from the 'test/testdat'
    directory
    """

    def _ignore_git(*args) -> str:  # pylint: disable=W0613
        """Copying the git directory is unnecessary so we ignore it

        This function is used by `shutil.copytree`"""
        return '.git'

    if os.path.exists('sandbox'):
        shutil.rmtree('sandbox')

    shutil.copytree('testdat', 'sandbox', ignore=_ignore_git)
    yield
    if os.path.exists('sandbox'):
        shutil.rmtree('sandbox')


""" Type representing the Iterator returned by `os.walk`"""
Walk = Iterator[Tuple[str, List[str], List[str]]]


@pytest.fixture(scope='session', autouse=True)
def initial_test_state() -> Walk:
    """Stores the initial state of the test data directory"""
    return os.walk('sandbox')


@pytest.fixture
def example_file(
        dat: Dict[str, Dict[str, List[str]]],
        tmpdir: py.path.local)-> Callable[[str, str], str]:
    """Return an example file generating fixture function"""

    def _get_example_file(filetype: str, size: str)-> str:
        """Return an example file of the requested filetype and size

        Parameters
        ----------
        filetype
            Name of a type of bioinformatics file
        size : {'tiny', 'small'}
            The approximate size of the output file. Tiny files have ~10
            entries, small files have ~1000-10000.
        """
        if filetype == 'fasta':
            outfile = dat[size]['assemblies'][0]
        elif filetype in ['fastq', 'fwd_fastq']:
            outfile = dat[size]['fwd_reads'][0]
        elif filetype == 'rev_fastq':
            outfile = dat[size]['rev_reads'][0]
        elif filetype == 'fai':
            outfile = dat[size]['faindices'][0]
        elif filetype == 'sam':
            outfile = dat[size]['sam'][0]
        elif filetype == 'bam':
            outfile = dat[size]['bam'][0]
        elif filetype == 'gz':
            outfile = dat[size]['zipped_fwd_reads'][0]
        else:
            outfile = gen_tmp(empty=False, directory=tmpdir, suffix='.foo')
        return outfile

    return _get_example_file


"""Generic type variable for the Filetyped classes Biofile and RuffusTask"""
T = TypeVar('T', fmruffus.RuffusTask, biofile.Biofile)


@pytest.fixture
def instance_of(example_file: Callable[[str, str], str]):
    """Given a class name, return an instance of the task

    Only works for classes with the class attributes input_type and/or
    output_type. Extra parameters cannot be passed. Designed for testing groups
    of closely related classes which are all initialized using the same basic
    process but with different types of input files.
    """
    def _make_test_instance(
            cls: Type[T],
            size: str)-> T:

        input_example = [example_file(t, size) for t in cls.input_type]
        input_example = fmlist.flatten(input_example)

        try:
            if len(input_example) == 1 and cls.output_type == ['']:
                output_example = [fmpaths.remove_suffix(input_example[0])]
            else:
                # If we have multiple inputs, the output suffix is added to
                # the first input as in ruffus
                input_prefix = fmpaths.remove_suffix(input_example[0])
                output_example = []
                for typ in cls.output_type:
                    output_example.append(
                        fmpaths.add_suffix(input_prefix, '.' + typ))

            if len(input_example) == 1:
                input_example = input_example[0]  # type: ignore
            if len(output_example) == 1:
                output_example = output_example[0]  # type: ignore
            return cls(input_example, output_example)  # type: ignore
        except (AttributeError, TypeError):
            return cls(*input_example)   # type: ignore
    return _make_test_instance


@pytest.fixture()
def nested_dir(tmpdir)-> str:
    """Create a set of nested directories and files inside a temp directory

    Returns
    -------
    Path of the temporary directory.
    """

    subdirs = ['foo', 'bar', 'car']
    subdirs = [os.path.join(str(tmpdir), d) for d in subdirs]
    for d in subdirs:
        os.makedirs(d)
        subdir_contents = ['a.x', 'b.y']
        subdir_contents = [os.path.join(d, f) for f in subdir_contents]
        for f in subdir_contents:
            open(f, 'a').close()
    return os.path.abspath(str(tmpdir))


@pytest.fixture(scope='session')
def dat() -> Dict[str, Dict[str, List[str]]]:
    """Create a dictionary of test data

    Assumes test directory is structured such that all test data is stored in
    the test/sandbox directory. test/sandbox can contain any number of
    directories which can each also contain any number of subdirectories. Each
    of these subdirectories contains a few test data files of the same type or
    use case. `dat` represents this structure as a dictionary of dictionaries.

    Returns
    -------
    A two level nested dictionary
    """
    # List contents of top level directories
    subdirs = fmpaths.listdirs('sandbox')
    dat = {}
    for d in subdirs:
        base = os.path.basename(d)
        dat[base] = fmpaths.contents_to_dict(d)
    return dat
