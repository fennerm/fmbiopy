"""Set of functions to aid in testing

Pytest functions must be imported explicitely.
"""
from os import chdir
from collections import namedtuple
from shutil import (
    copytree,
    rmtree,
)
from tempfile import NamedTemporaryFile
from uuid import uuid4

from plumbum import (
    FG,
    local,
)
from plumbum.cmd import (
    git,
    picard,
    samtools,
)
from pytest import fixture

from fmbiopy.fmpaths import (
    as_dict,
    as_paths,
    create_all,
    is_empty,
    listdirs,
    remove_all,
    root,
    silent_remove,
)


@fixture
def absolute_nonexist_paths(tmpdir, relative_nonexist_paths):
    """Generate a list of absolute paths to nonexistant paths"""
    return [tmpdir / p.name for p in relative_nonexist_paths]


@fixture
def absolute_exist_paths(tmpdir, randstrs):
    """Generate a list of absolute paths which exist"""
    paths = [tmpdir / x for x in randstrs(3)]
    create_all(paths)
    yield paths
    remove_all(paths, silent=True)


@fixture
def absolute_some_exist_paths(absolute_exist_paths, absolute_nonexist_paths):
    """Generate a list of paths, half of which exist"""
    return absolute_exist_paths + absolute_nonexist_paths


def assert_script_produces_files(script, args, output, redirect=None,
                                 empty_ok=False, outdir=None):
    """Assert that a script with given command line args produces expected files

    Parameters
    ----------
    script : str
        Path to the script
    args : List[str]
        List of command line arguments
    output: List[str] or List[plumbum.LocalPath]
        List of output files
    redirect: str or plumbum.LocalPath, optional
        If defined, redirect the stdout of the script to the given file.
    empty_ok : bool
        If True, output files are valid even if they are empty
    outdir: str or plumbum.LocalPath, optional
        If given, the output filenames are relative to this directory
    """
    execute = local[script]
    command = execute.__getitem__(args)

    if redirect:
        (command > redirect)()
    else:
        command()

    for f in output:
        if outdir:
            f = local.path(outdir) / f
        else:
            f = local.path(f)
        assert f.exists()
        if not empty_ok:
            assert not is_empty(f)


@fixture
def bowtie2_suffixes():
    """A list of the suffixes added by bowtie2-build"""
    return list(
        ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2'])


@fixture
def cd(tmpdir, startdir):
    """Change directory before running test"""
    chdir(str(tmpdir))
    yield
    chdir(str(startdir))


@fixture(scope='session')
def dat(sandbox):
    """Create a dictionary of test data

    Assumes test directory is structured such that all test data is stored in
    the test/sandbox directory. test/sandbox can contain any number of
    directories which can each also contain any number of subdirectories. Each
    of these subdirectories contains a few test data files of the same type or
    use case. `dat` represents this structure as a dictionary of dictionaries.

    Returns
    -------
    Dict[plumbum.LocalPath]
        A two level nested dictionary
    """
    # List contents of top level directories
    subdirs = listdirs(sandbox)
    dat = {}
    for d in subdirs:
        dat[d.name] = as_dict(d)
    return dat


@fixture
def double_suffixed_path(gen_tmp, tmpdir):
    """Generate a path with a two part suffix"""
    path = gen_tmp(empty=False, directory=tmpdir, suffix='.foo.bar')
    yield path
    silent_remove(path)


@fixture
def empty_list():
    """An empty list"""
    return []


@fixture
def empty_path(gen_tmp, tmpdir):
    """Generate an empty but existing path"""
    path = gen_tmp(directory=tmpdir, empty=True)
    yield path
    silent_remove(path)


@fixture()
def full_dir(tmpdir):
    """Create a temporary directory with some misc. temporary files"""
    paths = [tmpdir / name for name in ['a.x', 'b.y', 'b.x']]
    create_all(paths)
    yield tmpdir
    remove_all(paths, silent=True)


@fixture(scope='session')
def gen_tmp(sandbox):
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
    plumbum.LocalPath
    The path to the created temporary file
    """

    def _gen_tmp(empty=True, suffix='', directory=sandbox):

        tmpfile = local.path(
            NamedTemporaryFile(
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
    return local.path('.'.join([str(randpath()), randstr()[0:3], 'gz']))


@fixture(name='home')
def gen_home():
    """Get the path to the home directory"""
    return local.env.home()


@fixture(scope='session', autouse=True)
def load_sandbox(update_testdat, sandbox, testdat):
    """Copy all test data files to the sandbox for the testing session

    Yields
    ------
    A temporary directory with the contents copied from the 'test/testdat'
    directory
    """

    def _ignore_git(*args):
        """Copying the git directory is unnecessary so we ignore it

        This closure is used by `shutil.copytree`"""
        return '.git'

    if sandbox.exists():
        rmtree(str(sandbox), ignore_errors=True)

    copytree(str(testdat), str(sandbox), ignore=_ignore_git)
    (sandbox / '__init__.py').touch()
    yield sandbox
    if sandbox.exists():
        rmtree(str(sandbox), ignore_errors=True)


@fixture
def mixed_absolute_relative_paths(absolute_nonexist_paths,
                                  relative_nonexist_paths):
    """Generate a list of paths, half of which are absolute, half are not"""
    return absolute_nonexist_paths + relative_nonexist_paths


@fixture()
def nested_dir(tmpdir):
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
        rmtree(d, ignore_errors=True)


@fixture
def nonempty_path(gen_tmp, tmpdir):
    """Generate a nonempty path"""
    path = gen_tmp(empty=False, directory=tmpdir)
    yield path
    silent_remove(path)


@fixture
def nonempty_paths(gen_tmp, tmpdir):
    """Generate a list of nonempty paths"""
    paths = [gen_tmp(empty=False, directory=tmpdir) for i in range(3)]
    yield paths
    remove_all(paths, silent=True)


@fixture
def nonexistant_parent(randstr, tmpdir):
    """Generate a path for which the parent doesn't exist"""
    path = tmpdir / randstr() / randstr()
    yield path
    silent_remove(path)


@fixture()
def nonexistant_path(randstr, tmpdir):
    """Generate a nonexistant path"""
    path = tmpdir / randstr()
    yield path
    silent_remove(path)


@fixture(params=[
    'empty_path', 'nonexistant_path', 'nonempty_path', 'nonexistant_parent',
    'double_suffixed_path', 'symlink', 'tmpdir'
])
def poss_paths(request, empty_path, nonexistant_path, nonexistant_parent,
               nonempty_path, suffixed_path, double_suffixed_path, symlink,
               tmpdir):
    """Generate various kinds of possible valid `plumbum.LocalPath`s

    Returns
    -------
    A tuple of the form (name, value) where name is the name of the fixture
    """
    return (request.param, eval(request.param))


@fixture(params=[
    'absolute_exist_paths', 'absolute_nonexist_paths', 'empty_list',
    'relative_nonexist_paths', 'absolute_some_exist_paths',
    'mixed_absolute_relative_paths', 'nonempty_paths'
])
def poss_path_lists(request, absolute_exist_paths, absolute_nonexist_paths,
                    absolute_some_exist_paths, empty_list,
                    mixed_absolute_relative_paths, relative_nonexist_paths,
                    nonempty_paths):
    """Generate various kinds of possible valid lists of `LocalPath`s

    Returns
    -------
    A tuple of the form (name, value) where type is the name of the fixture
    """
    return (request.param, eval(request.param))


@fixture
def randpath(randstr, tmpdir):
    """Return a randomly generated nonexistant path"""

    def _gen_randpath():
        return tmpdir / randstr()

    return _gen_randpath


@fixture(scope='session')
def randstr():
    """Generate a unique random string"""

    def _get_rand_str():
        return uuid4().hex.replace('.', '')

    return _get_rand_str


@fixture(scope='session')
def randsuffix(randstr):
    """Generate a unique random suffix"""

    def get_rand_suffix():
        """test"""
        return '.' + randstr()

    return get_rand_suffix


@fixture(scope='session')
def randstrs(randstr):
    """Generate n unique random strings"""

    def _get_randstrs(n):
        return [randstr() for i in range(n)]

    return _get_randstrs


@fixture
def relative_nonexist_paths(randstrs, tmpdir):
    """Generate a list of relative paths"""
    paths = [local.path(tmpdir.name) / x for x in randstrs(3)]
    return paths


@fixture(scope='session')
def sandbox(testdir):
    """Path to the sandbox directory"""
    return testdir / 'sandbox'


@fixture(scope='session')
def small(sandbox):
    """Path to the 'small' subdirectory of sandbox"""
    return sandbox / 'small'


@fixture(scope='session')
def startdir():
    """The directory from which testing was started"""
    return local.cwd


@fixture
def suffixed_path(gen_tmp, tmpdir):
    """Generate a nonempty path with a suffix"""
    path = gen_tmp(empty=False, directory=tmpdir, suffix='.foo')
    yield path
    silent_remove(path)


@fixture
def suffixed_paths(gen_tmp, tmpdir, randsuffix):
    """Generate a list of nonempty suffixed paths"""
    suffixes = [randsuffix() for i in range(3)]
    paths = [
        gen_tmp(empty=False, directory=tmpdir, suffix=suffixes[i])
        for i in range(3)
    ]
    tup = namedtuple('suffixed_paths', ['paths', 'suffixes'])
    tup.paths = paths
    tup.suffixes = suffixes
    yield tup
    remove_all(paths)


@fixture
def symlink(gen_tmp, randpath, tmpdir):
    """Produce a symlink"""
    target = gen_tmp(empty=False, directory=tmpdir)
    path = randpath()
    path.symlink_to(target)
    yield path
    remove_all([target, path], silent=True)


@fixture(scope='session')
def testdat(testdir):
    """Path to the testdat directory"""
    return testdir / 'testdat'


@fixture(scope='session')
def testdat_repo():
    """The SSH address of the testdat github repo"""
    return 'git@github.com:fennerm/testdat'


@fixture(scope='session')
def testdir():
    """LocalPath to the test directory"""
    possible = as_paths(['test', 'tests', '.'])
    for poss in possible:
        if poss.exists():
            return poss
    raise OSError('Unsupported test directory structure')


@fixture(scope='session')
def tiny(sandbox):
    """Path to the 'small' subdirectory of sandbox"""
    return sandbox / 'tiny'


@fixture(name="tiny_indexed_bam")
def gen_tiny_indexed_bam(dat):
    bam = dat["tiny"]["bam"][0]
    samtools("index", bam)
    index = local.path(bam + ".bai")
    yield bam
    index.delete()


@fixture
def tmpdir(sandbox, randstr):
    """`plumbum.LocalPath` version of `pytest` fixture"""
    path = sandbox / randstr()
    path.mkdir()
    yield path
    rmtree(path)


@fixture()
def unique_dir(sandbox):
    """Create a directory in sandbox with a unique name"""
    path = sandbox / randstr()
    path.mkdir()
    yield path
    rmtree(str(path))


@fixture(scope='session', autouse=True)
def update_testdat(testdir, testdat, testdat_repo):
    """Make sure that testdat is up to date"""

    if not testdat.exists():
        with local.cwd(testdir):
            git['clone', testdat_repo]()
    else:
        with local.cwd(testdat):
            git['pull']()


def validate_bam_file(bam_or_sam):
    picard("ValidateSamFile", "I=" + bam_or_sam, "MODE=SUMMARY",
           "IGNORE_WARNINGS=true", "iGNORE=MISSING_READ_GROUP")
