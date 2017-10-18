"""Set of functions to aid in testing

Pytest functions must be imported explicitely.
"""

from pathlib import Path
from shutil import (
        copytree,
        rmtree,
        )
from tempfile import NamedTemporaryFile
from typing import (
        Callable,
        Dict,
        Iterator,
        List,
        )
from uuid import uuid4

from py import path
from pytest import fixture

from fmbiopy.fmpaths import (
        as_dict,
        as_paths,
        listdirs,
        )

@fixture(scope='session')
def gen_tmp(sandbox: Path)-> Callable[[bool, str, Path], Path]:
    """Generate a named temporary file.
    imported = importlib.import_module(module, package)

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
    def _gen_tmp(
            empty: bool = True,
            suffix: str = '',
            directory: Path = sandbox)-> Path:

        tmpfile = Path(NamedTemporaryFile(
            delete=False, dir=str(directory), suffix=suffix).name)

        if not empty:
            with tmpfile.open('w') as f:
                f.write('foo')
        else:
            with tmpfile.open('w') as f:
                f.write('')
        return tmpfile

    return _gen_tmp


@fixture(scope='session')
def sandbox(testdir: Path)-> Path:
    """Path to the sandbox directory"""
    return testdir / 'sandbox'

@fixture(scope='session')
def small(sandbox: Path)-> Path:
    """Path to the 'small' subdirectory of sandbox"""
    return sandbox / 'small'

@fixture(scope='session')
def tiny(sandbox: Path)-> Path:
    """Path to the 'small' subdirectory of sandbox"""
    return sandbox / 'tiny'


@fixture(scope='session', autouse=True)
def load_sandbox(sandbox: Path, testdat: Path) -> Iterator[Path]:
    """Copy all test data files to the sandbox for the testing session

    Yields
    ------
    A temporary directory with the contents copied from the 'test/testdat'
    directory
    """

    def _ignore_git(*args) -> str:  # pylint: disable=W0613
        """Copying the git directory is unnecessary so we ignore it

        This closure is used by `shutil.copytree`"""
        return '.git'

    if sandbox.exists():
        rmtree(str(sandbox))

    copytree(str(testdat), str(sandbox), ignore=_ignore_git)
    sandbox.joinpath('__init__.py').touch()
    yield sandbox
    if sandbox.exists():
        rmtree(str(sandbox))

@fixture()
def unique_dir(sandbox):
    """Create a directory in sandbox with a unique name"""
    path = sandbox.joinpath(uuid4().hex)
    path.mkdir()
    yield path
    rmtree(str(path))


@fixture(scope='session')
def testdir()-> Path:
    """Path to the test directory"""
    possible = as_paths(['test', 'tests', '.'])
    for poss in possible:
        if poss.exists():
            return poss
    raise OSError('Unsupported test directory structure')



@fixture(scope='session')
def testdat(testdir)-> Path:
    """Path to the testdat directory"""
    return testdir / 'testdat'


@fixture(scope='session')
def example_file(
        gen_tmp: Callable[[bool, str, Path], Path],
        dat: Dict[str, Dict[str, List[str]]])-> Callable[[str, str], Path]:
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
            return gen_tmp(empty=False, suffix='.foo')
        return Path(outfile)

    return _get_example_file


@fixture(autouse=True)
def tmpdir(tmpdir: path.local)-> Path:
    """`pathlib.Path` version of `pytest` fixture"""
    return Path(str(tmpdir))


@fixture()
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


@fixture()
def full_dir(tmpdir: path.local)-> Path:
    """Create a temporary directory with some misc. temporary files"""
    paths = [tmpdir / name for name in ['a.x', 'b.y', 'b.x']]
    for p in paths:
        p.touch()
    return tmpdir


@fixture(scope='session')
def dat(sandbox) -> Dict[str, Dict[str, List[Path]]]:
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
    subdirs = listdirs(sandbox)
    dat = {}
    for d in subdirs:
        dat[d.name] = as_dict(d)
    return dat
