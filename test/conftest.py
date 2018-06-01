"""pytest fixtures shared among tests"""
from __future__ import print_function
from os import chdir
from collections import namedtuple
from shutil import copytree, rmtree
from tempfile import NamedTemporaryFile
from uuid import uuid4

from Bio import SeqIO
from numpy.random import binomial
import pandas as pd
from plumbum import local
from plumbum.cmd import git, sambamba, samtools
from pytest import fixture

from fmbiopy.fmbio import (
    align_and_sort,
    count_reads,
    index_fasta,
    simulate_fasta,
)
from fmbiopy.fmpaths import (
    as_dict,
    is_empty,
    listdirs,
    remove_all,
    silent_remove,
)
from fmbiopy.fmsystem import capture_stdout
from test.helpers import trim


@fixture
def absolute_nonexist_paths(tmpdir, relative_nonexist_paths):
    """Generate a list of absolute paths to nonexistant paths"""
    return [tmpdir / p.name for p in relative_nonexist_paths]


@fixture
def absolute_exist_paths(tmpdir, randstrs):
    """Generate a list of absolute paths which exist"""
    paths = [tmpdir / x for x in randstrs(3)]
    for path in paths:
        path.touch()
    yield paths
    remove_all(paths, silent=True)


@fixture
def absolute_some_exist_paths(absolute_exist_paths, absolute_nonexist_paths):
    """Generate a list of paths, half of which exist"""
    return absolute_exist_paths + absolute_nonexist_paths


@fixture
def bowtie2_suffixes():
    """A list of the suffixes added by bowtie2-build"""
    return list(
        [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
    )


@fixture
def cd(tmpdir, startdir):
    """Change directory before running test"""
    chdir(str(tmpdir))
    yield
    chdir(str(startdir))


@fixture(scope="session")
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
    path = gen_tmp(empty=False, directory=tmpdir, suffix=".foo.bar")
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
    paths = [tmpdir / name for name in ["a.x", "b.y", "b.x"]]
    for path in paths:
        path.touch()
    yield tmpdir
    remove_all(paths, silent=True)


@fixture(scope="session")
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

    def _gen_tmp(empty=True, suffix="", directory=sandbox):

        tmpfile = local.path(
            NamedTemporaryFile(
                delete=False, dir=str(directory), suffix=suffix
            ).name
        )

        if not empty:
            with tmpfile.open("w") as f:
                f.write("foo")
        else:
            with tmpfile.open("w") as f:
                f.write("")
        return tmpfile

    return _gen_tmp


@fixture
def gzipped_path(randpath, randstr):
    """Return a gzipped file name"""
    return local.path(".".join([str(randpath()), randstr()[0:3], "gz"]))


@fixture
def home():
    """Get the path to the home directory"""
    return local.env.home()


@fixture(scope="session", autouse=True)
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
        return ".git"

    if sandbox.exists():
        rmtree(str(sandbox), ignore_errors=True)

    copytree(str(testdat), str(sandbox), ignore=_ignore_git)
    (sandbox / "__init__.py").touch()
    yield sandbox
    if sandbox.exists():
        rmtree(str(sandbox), ignore_errors=True)


@fixture
def mixed_absolute_relative_paths(
    absolute_nonexist_paths, relative_nonexist_paths
):
    """Generate a list of paths, half of which are absolute, half are not"""
    return absolute_nonexist_paths + relative_nonexist_paths


@fixture
def nested_dir(tmpdir):
    """Create a set of nested directories and files inside a temp directory

    Yields
    ------
    LocalPath
    Path of the temporary directory.
    """

    subdir_names = ["foo", "bar", "car"]
    subdirs = [tmpdir / d for d in subdir_names]
    for d in subdirs:
        d.mkdir()
        subdir_contents = [d / f for f in ["a.x", "b.y"]]
        for f in subdir_contents:
            f.open("a").close()
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


@fixture
def nonexistant_path(randstr, tmpdir):
    """Generate a nonexistant path"""
    path = tmpdir / randstr()
    yield path
    silent_remove(path)


@fixture(
    params=[
        "empty_path",
        "nonexistant_path",
        "nonempty_path",
        "nonexistant_parent",
        "double_suffixed_path",
        "symlink",
        "tmpdir",
    ]
)
def poss_paths(
    request,
    empty_path,
    nonexistant_path,
    nonexistant_parent,
    nonempty_path,
    suffixed_path,
    double_suffixed_path,
    symlink,
    tmpdir,
):
    """Generate various kinds of possible valid `plumbum.LocalPath`s

    Returns
    -------
    A tuple of the form (name, value) where name is the name of the fixture
    """
    return (request.param, eval(request.param))


@fixture(
    params=[
        "absolute_exist_paths",
        "absolute_nonexist_paths",
        "empty_list",
        "relative_nonexist_paths",
        "absolute_some_exist_paths",
        "mixed_absolute_relative_paths",
        "nonempty_paths",
    ]
)
def poss_path_lists(
    request,
    absolute_exist_paths,
    absolute_nonexist_paths,
    absolute_some_exist_paths,
    empty_list,
    mixed_absolute_relative_paths,
    relative_nonexist_paths,
    nonempty_paths,
):
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


@fixture(scope="session")
def randstr():
    """Generate a unique random string"""

    def _get_rand_str():
        return uuid4().hex.replace(".", "")

    return _get_rand_str


@fixture(scope="session")
def randsuffix(randstr):
    """Generate a unique random suffix"""

    def get_rand_suffix():
        """test"""
        return "." + randstr()

    return get_rand_suffix


@fixture(scope="session")
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


@fixture(scope="session")
def sandbox(testdir):
    """Path to the sandbox directory"""
    return testdir / "sandbox"


@fixture(scope="session")
def small(sandbox):
    """Path to the 'small' subdirectory of sandbox"""
    return sandbox / "small"


@fixture(scope="session")
def startdir():
    """The directory from which testing was started"""
    return local.cwd


@fixture
def suffixed_path(gen_tmp, tmpdir):
    """Generate a nonempty path with a suffix"""
    path = gen_tmp(empty=False, directory=tmpdir, suffix=".foo")
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
    tup = namedtuple("suffixed_paths", ["paths", "suffixes"])
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


@fixture(scope="session")
def testdat(testdir):
    """Path to the testdat directory"""
    return testdir / "testdat"


@fixture(scope="session")
def testdat_repo():
    """The SSH address of the testdat github repo"""
    return "git@github.com:fennerm/testdat"


@fixture(scope="session")
def testdir():
    """LocalPath to the test directory"""
    possible = [local.path(x) for x in ["test", "tests", "."]]
    for poss in possible:
        if poss.exists():
            return poss
    raise OSError("Unsupported test directory structure")


@fixture(scope="session")
def tiny(sandbox):
    """Path to the 'small' subdirectory of sandbox"""
    return sandbox / "tiny"


@fixture(scope="session")
def fasta(sandbox):
    output = {}
    output["fasta"] = sandbox / (uuid4().hex + ".fa")
    output["bt2"] = output["fasta"].with_suffix("")
    simulate_fasta(20, 500, output["fasta"])
    index_fasta(output["fasta"], "all")
    return output


@fixture(scope="session")
def nonindexed_fasta(sandbox, fasta):
    output_file = sandbox / (uuid4().hex + ".fa")
    fasta["fasta"].copy(output_file)
    return output_file


@fixture(scope="session")
def simulated_reads(sandbox, fasta):
    python2 = local["python2"]
    gen_reads = python2["test/lib/neat-genreads/genReads.py"]
    output_prefix = sandbox / uuid4().hex
    output = {}
    output["fwd"] = local.path(output_prefix + "_read1.fq")
    output["rev"] = local.path(output_prefix + "_read2.fq")
    output["bam"] = local.path(output_prefix + "_golden.bam")

    gen_reads[
        "-r",
        fasta["fasta"],
        "-R",
        "101",
        "-o",
        output_prefix,
        "--bam",
        "--pe",
        "300",
        "30",
    ]()
    return output


@fixture(scope="session")
def paired_trimmed_fastq(sandbox, simulated_reads):
    """Produce a test dataset with trimmed paired fastq files

    Returns a dict with 'fwd', 'rev' and 'unpaired' entries.
    """
    prefix = sandbox / uuid4().hex
    prob_removal = 0.05
    prob_trim = 0.5
    trim_interval = [1, 30]

    trimmed = {}
    trimmed["fwd"] = local.path(prefix + ".R1.fastq")
    trimmed["rev"] = local.path(prefix + ".R2.fastq")
    trimmed["unpaired"] = local.path(prefix + ".unpaired.fastq")

    with trimmed["fwd"].open("w") as out_fwd, trimmed["rev"].open(
        "w"
    ) as out_rev, trimmed["unpaired"].open("w") as out_unp:
        for reads in zip(
            SeqIO.parse(simulated_reads["fwd"], "fastq"),
            SeqIO.parse(simulated_reads["rev"], "fastq"),
        ):
            # Trim the reads
            reads = [trim(read, prob_trim, trim_interval) for read in reads]

            # Remove some of the reads
            is_removed = binomial(1, prob_removal, 2)
            if is_removed[0] and not is_removed[1]:
                out_unp.write(reads[1].format("fastq"))
            elif is_removed[1] and not is_removed[0]:
                out_unp.write(reads[0].format("fastq"))
            else:
                out_fwd.write(reads[0].format("fastq"))
                out_rev.write(reads[1].format("fastq"))
    for fastq in list(trimmed.values()):
        assert count_reads(fastq) > 0

    return trimmed


@fixture(scope="session")
def partial_fasta(sandbox, fasta):
    output = {}
    output["fasta"] = sandbox / (uuid4().hex + ".fa")
    output["bt2"] = output["fasta"].with_suffix("")
    with output["fasta"].open("w") as f:
        for i, record in enumerate(SeqIO.parse(fasta["fasta"], "fasta")):
            if i < 4:
                # Full sequences
                print(">" + record.id, file=f)
                print(record.seq, file=f)
            elif i > 3 and i < 7:
                # Trimmed sequences
                print(">" + record.id, file=f)
                print(record.seq[0:250], file=f)
    index_fasta(output["fasta"], "all")
    return output


@fixture(scope="session")
def mini_bams(sandbox, trimmed_bam):
    list_csomes = local["bin/list_csomes"]
    csomes = capture_stdout(list_csomes[trimmed_bam])
    output_bams = []
    for csome in csomes:
        output_bam = sandbox / (uuid4().hex + ".bam")
        output_bams.append(output_bam)
        (samtools["view", "-bh", trimmed_bam, csome] > output_bam)()
        assert count_reads(output_bam) > 0
    return output_bams


@fixture(scope="session")
def empty_bam(sandbox, untrimmed_bam):
    output_bam = sandbox / (uuid4().hex + ".bam")
    (sambamba["view", "-f", "bam", "-F", "ref_name==foo"] > output_bam)()
    assert count_reads(output_bam) == 0
    return output_bam


@fixture(scope="session")
def untrimmed_bam(sandbox, fasta, simulated_reads):
    output_bam = sandbox / (uuid4().hex + ".bam")
    align_and_sort(
        idx=fasta["bt2"],
        fastq1=simulated_reads["fwd"],
        fastq2=simulated_reads["rev"],
        preset="very-fast",
        output_bam=output_bam,
    )
    return output_bam


@fixture(scope="session")
def trimmed_bam(sandbox, partial_fasta, paired_trimmed_fastq):
    output_bam = sandbox / (uuid4().hex + ".bam")
    align_and_sort(
        idx=partial_fasta["bt2"],
        fastq1=paired_trimmed_fastq["fwd"],
        fastq2=paired_trimmed_fastq["rev"],
        preset="very-fast",
        output_bam=output_bam,
    )
    return output_bam


@fixture(scope="session")
def bam_with_orphans(sandbox, trimmed_bam):
    """Generate a .bam file with orphaned reads."""
    output_bam = sandbox / (uuid4().hex + ".bam")
    # Reads which match this criterion are filtered, leaving orphaned reads in
    # the bam
    filt = "not (first_of_pair and ref_id == 1 and position < 200)"
    (sambamba["view", "-f", "bam", "-F", filt, trimmed_bam] > output_bam)()
    return output_bam


@fixture
def indexed_bam_with_orphans(bam_with_orphans):
    """A bam file with orphan reads and a .bai index."""
    samtools("index", bam_with_orphans)
    return bam_with_orphans


@fixture(scope="session")
def trimmed_sam(sandbox, trimmed_bam):
    bam_to_sam = local["bin/bam_to_sam"]
    sam = sandbox / (uuid4().hex + ".sam")
    (bam_to_sam[trimmed_bam] > sam)()
    return sam


@fixture(scope="session")
def untrimmed_sam(sandbox, untrimmed_bam):
    bam_to_sam = local["bin/bam_to_sam"]
    sam = sandbox / (uuid4().hex + ".sam")
    (bam_to_sam[untrimmed_bam] > sam)()
    return sam


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


@fixture(scope="session", autouse=True)
def update_testdat(testdir, testdat, testdat_repo):
    """Make sure that testdat is up to date"""

    if not testdat.exists():
        with local.cwd(testdir):
            git["clone", testdat_repo]()
    else:
        with local.cwd(testdat):
            git["pull"]()


@fixture
def dataframe():
    """A dataframe with following structure:

       A  B
    1  0  0
    2  0  0

    """
    x = pd.DataFrame(index=[1, 2], columns=["A", "B"])
    x = x.fillna(0)
    return x


@fixture
def dataframe_header():
    return ["# foo\n", "\n"]
