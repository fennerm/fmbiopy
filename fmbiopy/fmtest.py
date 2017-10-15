"""Set of functions to aid in testing

Pytest functions must be imported explicitely.
"""

from pathlib import Path
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

import fmbiopy.fmpaths as fmpaths


def gen_tmp(
        empty: bool = True,
        suffix: str = '',
        directory: str = 'sandbox') -> Path:
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

    tmpfile = Path(tempfile.NamedTemporaryFile(
        delete=False, dir=directory, suffix=suffix).name)

    if not empty:
        with tmpfile.open('w') as f:
            f.write('foo')
    else:
        with tmpfile.open('w') as f:
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

    sandbox = Path('sandbox')
    if sandbox.exists():
        shutil.rmtree(sandbox.name)

    shutil.copytree('testdat', 'sandbox', ignore=_ignore_git)
    yield sandbox
    if sandbox.exists():
        shutil.rmtree(sandbox.name)


@pytest.fixture
def example_file(
        dat: Dict[str, Dict[str, List[str]]],
        tmpdir: py.path.local)-> Callable[[str, str], Path]:
    """Return an example file generating fixture function"""

    def _get_example_file(filetype: str, size: str)-> Path:
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
            return gen_tmp(empty=False, directory=tmpdir, suffix='.foo')
        return Path(outfile)

    return _get_example_file


@pytest.fixture(autouse=True)
def tmpdir(tmpdir: py.path.local)-> Path:
    return Path(tmpdir)


@pytest.fixture()
def nested_dir(tmpdir: Path)-> Path:
    """Create a set of nested directories and files inside a temp directory

    Returns
    -------
    Path of the temporary directory.
    """

    subdir_names = ['foo', 'bar', 'car']
    subdirs = [tmpdir / d for d in subdir_names]
    for d in subdirs:
        d.mkdir()
        subdir_contents = [d / f for f in ['a.x', 'b.y']]
        for f in subdir_contents:
            f.open('a').close()
    return Path(tmpdir)


@pytest.fixture()
def full_dir(tmpdir: py.path.local)-> Path:
    """Create a temporary directory with some misc. temporary files"""
    paths = [tmpdir / name for name in ['a.x', 'b.y', 'b.x']]
    for path in paths:
        path.touch()
    return tmpdir


@pytest.fixture(scope='session')
def dat() -> Dict[str, Dict[str, List[Path]]]:
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
    subdirs = fmpaths.listdirs(Path('sandbox'))
    dat = {}
    for d in subdirs:
        dat[d.name] = fmpaths.as_dict(d)
    return dat
