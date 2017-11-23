"""Set of functions to aid in testing

Pytest functions must be imported explicitely.
"""
from os import chdir
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
        NamedTuple,
        Tuple,
        Type,
        )
from uuid import uuid4

import biofile as bf
from fmbiopy.fmpaths import (
        as_dict,
        as_paths,
        create_all,
        listdirs,
        remove_all,
        root,
        silent_remove,
        )
from plumbum.cmd import git
from plumbum import (
        local,
        LocalPath,
        )
from pytest import fixture
from _pytest.fixtures import SubRequest


"""The type of the gen_tmp fixture"""
GenTmpType = Callable[[bool, str, LocalPath], LocalPath]


@fixture
def absolute_nonexist_paths(
        tmpdir: LocalPath,
        relative_nonexist_paths: List[LocalPath])-> List[LocalPath]:
    """Generate a list of absolute paths to nonexistant paths"""
    return [tmpdir / p.name for p in relative_nonexist_paths]


@fixture
def absolute_exist_paths(
        tmpdir: LocalPath,
        randstrs: Callable[[int], List[str]])-> Iterator[List[LocalPath]]:
    """Generate a list of absolute paths which exist"""
    paths = [tmpdir / x for x in randstrs(3)]
    create_all(paths)
    yield paths
    remove_all(paths, silent=True)


@fixture
def absolute_some_exist_paths(
        absolute_exist_paths: List[LocalPath],
        absolute_nonexist_paths: List[LocalPath])-> List[LocalPath]:
    """Generate a list of paths, half of which exist"""
    return absolute_exist_paths + absolute_nonexist_paths

@fixture
def bowtie2_suffixes():
    """A list of the suffixes added by bowtie2-build"""
    return list(['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
                   '.rev.2.bt2'])

@fixture
def cd(tmpdir, startdir)-> Iterator[None]:
    """Change directory before running test"""
    chdir(str(tmpdir))
    yield
    chdir(str(startdir))


@fixture(scope='session')
def dat(sandbox) -> Dict[str, Dict[str, List[LocalPath]]]:
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



@fixture
def double_suffixed_path(
        gen_tmp: GenTmpType,
        tmpdir: LocalPath)-> Iterator[LocalPath]:
    """Generate a path with a two part suffix"""
    path = gen_tmp(empty=False, directory=tmpdir, suffix='.foo.bar')
    yield path
    silent_remove(path)

@fixture
def empty_list()-> List:
    """An empty list"""
    return []


@fixture
def empty_path(gen_tmp: GenTmpType, tmpdir: LocalPath)-> Iterator[LocalPath]:
    """Generate an empty but existing path"""
    path = gen_tmp(directory=tmpdir, empty=True)
    yield path
    silent_remove(path)


@fixture(scope='session')
def example_file(
        gen_tmp: GenTmpType,
        dat: Dict[str, Dict[str, List[LocalPath]]]
        )-> Callable[[bf.Biofile, str], LocalPath]:
    """Return an example file generating fixture function"""

    def _get_example_file(filetype: Type[bf.Biofile], size: str)-> LocalPath:
        """Return an example file of the requested filetype and size

        Parameters
        ----------
        filetype
            Name of a type of bioinformatics file
        size : {'tiny', 'small'}
            The approximate size of the output file. Tiny files have ~10
            entries, small files have ~1000-10000.

        Returns
        -------
        LocalPath
            A path to an example file of the requested type.
        """
        try:
            cls = filetype
        except AttributeError:
            cls = bf.type_to_class(filetype)

        if cls.__name__ == 'Fasta':
            outfile = dat[size]['assemblies'][0]
        elif cls.__name__ in ['Fastq', 'FwdFastq']:
            outfile = dat[size]['fwd_reads'][0]
        elif cls.__name__ == 'RevFastq':
            outfile = dat[size]['rev_reads'][0]
        elif cls.__name__ == 'SamtoolsFAIndex':
            outfile = dat[size]['faindices'][0]
        elif cls.__name__ == 'Sam':
            outfile = dat[size]['sam'][0]
        elif cls.__name__ == 'Bam':
            outfile = dat[size]['bam'][0]
        elif cls.__name__ == 'Adapters':
            outfile = dat[size]['adapters'][0]
        elif cls.__name__ == 'Gzipped':
            outfile = dat[size]['zipped_fwd_reads'][0]
        elif cls.__name__ == 'CentrifugeDB':
            outfile = root(dat[size]['centrifuge_idx'][0])
        elif cls.__name__ == 'CentrifugeOutput':
            outfile = dat[size]['centrifuge_output'][0]
        elif cls.__name__ == 'UnpairedFastq':
            outfile = dat[size]['fwd_reads'][1]
        elif cls.__name__ == 'Biofile':
            outfile = str(gen_tmp(empty=False, suffix='.foo'))
        else:
            outfile = str(gen_tmp(empty=False, suffix=cls.extensions[0]))
        return local.path(outfile)

    return _get_example_file


@fixture()
def full_dir(tmpdir: LocalPath)-> Iterator[LocalPath]:
    """Create a temporary directory with some misc. temporary files"""
    paths = [tmpdir / name for name in ['a.x', 'b.y', 'b.x']]
    create_all(paths)
    yield tmpdir
    remove_all(paths, silent=True)


@fixture(scope='session')
def gen_tmp(sandbox: LocalPath)-> GenTmpType:
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
            directory: LocalPath = sandbox)-> LocalPath:

        tmpfile = LocalPath(NamedTemporaryFile(
            delete=False, dir=str(directory), suffix=suffix).name)

        if not empty:
            with tmpfile.open('w') as f:
                f.write('foo')
        else:
            with tmpfile.open('w') as f:
                f.write('')
        return tmpfile

    return _gen_tmp


@fixture
def gzipped_path(randpath, randstr):
    """Return a gzipped file name"""
    return LocalPath('.'.join([str(randpath()), randstr()[0:3], 'gz']))

@fixture(name='home')
def gen_home():
    """Get the path to the home directory"""
    return LocalPath.home()

@fixture(scope='session', autouse=True)
def load_sandbox(
        update_testdat: None,
        sandbox: LocalPath,
        testdat: LocalPath) -> Iterator[LocalPath]:
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
    (sandbox / '__init__.py').touch()
    yield sandbox
    if sandbox.exists():
        rmtree(str(sandbox))


@fixture
def mixed_absolute_relative_paths(
        absolute_nonexist_paths: List[LocalPath],
        relative_nonexist_paths: List[LocalPath])-> List[LocalPath]:
    """Generate a list of paths, half of which are absolute, half are not"""
    return absolute_nonexist_paths + relative_nonexist_paths


@fixture()
def nested_dir(tmpdir: LocalPath)-> Iterator[LocalPath]:
    """Create a set of nested directories and files inside a temp directory

    Yields
    ------
    LocalPath
        Path of the temporary directory.
    """

    subdir_names = ['foo', 'bar', 'car']
    subdirs = [tmpdir / d for d in subdir_names]
    for d in subdirs:
        d.mkdir()
        subdir_contents = [d / f for f in ['a.x', 'b.y']]
        for f in subdir_contents:
            f.open('a').close()
    yield tmpdir
    for d in subdirs:
        rmtree(d)


@fixture
def nonempty_path(gen_tmp: GenTmpType, tmpdir: LocalPath)-> Iterator[LocalPath]:
    """Generate a nonempty path"""
    path = gen_tmp(empty=False, directory=tmpdir)
    yield path
    silent_remove(path)


@fixture
def nonempty_paths(
        gen_tmp: GenTmpType,
        tmpdir: LocalPath)-> Iterator[List[LocalPath]]:
    """Generate a list of nonempty paths"""
    paths = [gen_tmp(empty=False, directory=tmpdir) for i in range(3)]
    yield paths
    remove_all(paths, silent=True)


@fixture
def nonexistant_parent(
        randstr: Callable[[], str],
        tmpdir: LocalPath)-> Iterator[LocalPath]:
    """Generate a path for which the parent doesn't exist"""
    path = tmpdir / randstr() / randstr()
    yield path
    silent_remove(path)


@fixture()
def nonexistant_path(
        randstr: Callable[[], str],
        tmpdir: LocalPath)-> Iterator[LocalPath]:
    """Generate a nonexistant path"""
    path = tmpdir / randstr()
    yield path
    silent_remove(path)


@fixture(params=[
    'empty_path', 'nonexistant_path', 'nonempty_path', 'nonexistant_parent',
    'double_suffixed_path', 'symlink', 'tmpdir'])
def poss_paths(
        request: SubRequest,
        empty_path: LocalPath,
        nonexistant_path: LocalPath,
        nonexistant_parent: LocalPath,
        nonempty_path: LocalPath,
        suffixed_path: LocalPath,
        double_suffixed_path: LocalPath,
        symlink: LocalPath,
        tmpdir: LocalPath,
        )-> Tuple[str, LocalPath]:
    """Generate various kinds of possible valid `LocalPath`s

    Returns
    -------
    A tuple of the form (name, value) where name is the name of the fixture
    """
    return (request.param, eval(request.param))


@fixture(params=[
    'absolute_exist_paths', 'absolute_nonexist_paths', 'empty_list',
    'relative_nonexist_paths', 'absolute_some_exist_paths',
    'mixed_absolute_relative_paths', 'nonempty_paths'])
def poss_path_lists(
        request: SubRequest,
        absolute_exist_paths: List[LocalPath],
        absolute_nonexist_paths: List[LocalPath],
        absolute_some_exist_paths: List[LocalPath],
        empty_list: List,
        mixed_absolute_relative_paths: List[LocalPath],
        relative_nonexist_paths: List[LocalPath],
        nonempty_paths: List[LocalPath],
        )-> Tuple[str, LocalPath]:
    """Generate various kinds of possible valid lists of `LocalPath`s

    Returns
    -------
    A tuple of the form (name, value) where type is the name of the fixture
    """
    return (request.param, eval(request.param))


""" Type variable for `randpath`"""
RandPathType = Callable[[], LocalPath]


@fixture
def randpath(randstr: Callable[[], str], tmpdir: LocalPath)-> RandPathType:
    """Return a randomly generated nonexistant path"""
    def _gen_randpath()-> LocalPath:
        return tmpdir / randstr()
    return _gen_randpath


@fixture(scope='session')
def randstr()-> Callable[[], str]:
    """Generate a unique random string"""
    def _get_rand_str()-> str:
        return uuid4().hex.replace('.', '')
    return _get_rand_str

@fixture(scope='session')
def randsuffix(randstr)-> Callable[[], str]:
    """Generate a unique random suffix"""
    def get_rand_suffix()-> str:
        """test"""
        return '.' + randstr()
    return get_rand_suffix


@fixture(scope='session')
def randstrs(randstr: Callable[[], str])-> Callable[[int], List[str]]:
    """Generate n unique random strings"""
    def _get_randstrs(n: int)-> List[str]:
        return [randstr() for i in range(n)]
    return _get_randstrs


@fixture
def relative_nonexist_paths(
        randstrs: Callable[[int], List[str]],
        tmpdir: LocalPath)-> List[LocalPath]:
    """Generate a list of relative paths"""
    paths = [LocalPath(tmpdir.name) / x for x in randstrs(3)]
    return paths


@fixture(scope='session')
def sandbox(testdir: LocalPath)-> LocalPath:
    """Path to the sandbox directory"""
    return testdir / 'sandbox'


@fixture(scope='session')
def small(sandbox: LocalPath)-> LocalPath:
    """Path to the 'small' subdirectory of sandbox"""
    return sandbox / 'small'


@fixture(scope='session')
def startdir()-> LocalPath:
    """The directory from which testing was started"""
    return local.cwd

@fixture
def suffixed_path(gen_tmp: GenTmpType, tmpdir: LocalPath)-> Iterator[LocalPath]:
    """Generate a nonempty path with a suffix"""
    path = gen_tmp(empty=False, directory=tmpdir, suffix='.foo')
    yield path
    silent_remove(path)

@fixture
def suffixed_paths(
        gen_tmp: GenTmpType,
        tmpdir: LocalPath,
        randsuffix: Callable[[], str]
        )-> Iterator[NamedTuple]:
    """Generate a list of nonempty suffixed paths"""
    suffixes = [randsuffix() for i in range(3)]
    paths = [
            gen_tmp(empty=False, directory=tmpdir, suffix=suffixes[i])
            for i in range(3)]
    tup = NamedTuple('suffixed_paths', [
        ('paths', List[LocalPath]), ('suffixes', List[str])])
    tup.paths = paths
    tup.suffixes = suffixes
    yield tup
    remove_all(paths)


@fixture
def symlink(
        gen_tmp: GenTmpType,
        randpath: RandPathType,
        tmpdir: LocalPath,
        )-> Iterator[LocalPath]:
    """Produce a symlink"""
    target = gen_tmp(empty=False, directory=tmpdir)
    path = randpath()
    path.symlink_to(target)
    yield path
    remove_all([target, path], silent=True)


@fixture(scope='session')
def testdat(testdir)-> LocalPath:
    """Path to the testdat directory"""
    return testdir / 'testdat'

@fixture(scope='session')
def testdat_repo()-> str:
    """The SSH address of the testdat github repo"""
    return 'git@github.com:fennerm/testdat'


@fixture(scope='session')
def testdir()-> LocalPath:
    """LocalPath to the test directory"""
    possible = as_paths(['test', 'tests', '.'])
    for poss in possible:
        if poss.exists():
            return poss
    raise OSError('Unsupported test directory structure')


@fixture(scope='session')
def tiny(sandbox: LocalPath)-> LocalPath:
    """Path to the 'small' subdirectory of sandbox"""
    return sandbox / 'tiny'


@fixture
def tmpdir(
        sandbox: LocalPath,
        randstr: Callable[[], str])-> Iterator[LocalPath]:
    """`plumbum.LocalPath` version of `pytest` fixture"""
    path = sandbox / randstr()
    path.mkdir()
    yield path
    rmtree(path)


@fixture()
def unique_dir(sandbox: LocalPath)-> Iterator[LocalPath]:
    """Create a directory in sandbox with a unique name"""
    path = sandbox / randstr()
    path.mkdir()
    yield path
    rmtree(str(path))


@fixture(scope='session', autouse=True)
def update_testdat(
        testdir: LocalPath, testdat: LocalPath, testdat_repo: str)-> None:
    """Make sure that testdat is up to date"""

    if not testdat.exists():
        with local.cwd(testdir):
            git['clone', testdat_repo]()
    else:
        with local.cwd(testdat):
            git['pull']()
