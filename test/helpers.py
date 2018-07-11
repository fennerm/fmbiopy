"""Test helper functions."""
from functools import partial
import os
from random import randint
from uuid import uuid4

from pandas import DataFrame
from numpy.random import binomial
from plumbum import local
from plumbum.cmd import picard
import wrapt

from fmbiopy.obj import get_param_names, replace_param_sig
from fmbiopy.paths import is_empty

SANDBOX = local.path("test/sandbox")


def assert_script_produces_files(
    script, args, output, redirect=None, empty_ok=False, outdir=None
):
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


def trim(read, prob_trim, trim_interval):
    if binomial(1, prob_trim) == 1:
        trimmed_bases = randint(trim_interval[0], trim_interval[1])
        return read[:-trimmed_bases]
    return read


def validate_bam_file(bam_or_sam):
    picard(
        "ValidateSamFile",
        "I=" + bam_or_sam,
        "MODE=SUMMARY",
        "IGNORE_WARNINGS=true",
        "iGNORE=MISSING_READ_GROUP",
    )


def assert_df_equals(df, query):
    """Check if a DataFrame is equal to a DataFrame or dictionary."""
    if isinstance(query, dict):
        query = DataFrame.from_dict(query)
    assert df.to_dict() == query.to_dict()


def gen_reads(
    fasta, output_dir, bam_output=True, vcf_output=False, mutation_rate=0
):
    """Wrapper around NEAT gen-reads script."""
    output_prefix = output_dir / uuid4().hex
    output = {
        "prefix": output_prefix,
        "fwd": local.path(output_prefix + "_read1.fq"),
        "rev": local.path(output_prefix + "_read2.fq"),
        "bam": local.path(output_prefix + "_golden.bam"),
        "vcf": local.path(output_prefix + "_golden.vcf"),
    }
    gen_reads_bin = local["python2"]["test/lib/neat-genreads/genReads.py"]
    gen_reads_args = [
        "-r",
        fasta,
        "-R",
        "101",
        "-o",
        output_prefix,
        "-M",
        mutation_rate,
        "--pe",
        "300",
        "30",
    ]
    if bam_output:
        gen_reads_args.append("--bam")
    if vcf_output:
        gen_reads_args.append("--vcf")
    gen_reads_bin.__getitem__(gen_reads_args)()
    return output


def file_generator(
    wrapped=None,
    ids=["file"],
    names=[uuid4().hex],
    suffixes=[""],
    dirs=[SANDBOX],
    properties=None,
):
    """Decorator which automates setup and return for file generation functions.

    The decorator fulfills 3 tasks:
        1. Generating required temporary file names.
        2. Creating an output dictionary (`D`) which lists file locations.
        3. Modifies the return of the wrapped function to instead return `D`.

    Wrapped functions must contain a parameter named `meta`, this parameter will
    be replaced by `D` at runtime, allowing the function to access the generated
    filenames.

    Parameters
    ----------
    wrapped: Callable
        Should not be defined manually, will be passed automatically by python.
    ids: Iterable[str]
        A list of ids for output files. These will be the keys for the output
        dictionary.
    names: Iterable[str]
        A list of output file basenames.
    suffixes: Iterable[str]
        A list of file suffixes for output files.
    dirs: Iterable[PathLike]
        A list of output directories.
    properties: Dict
        Dictionary of extra metadata to be included in output dictionary.

    """
    # warning: this function gets pretty hairy

    if wrapped is None:
        # Happens when some optional params are defined
        # See wrapt docs for reasoning
        return partial(
            file_generator,
            ids=ids,
            names=names,
            suffixes=suffixes,
            dirs=dirs,
            properties=properties,
        )

    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        # Wrapper function with decorator arguments included implicitely as
        # variables in outer scope.

        # First generate the metadata which will be passed into the decorated
        # function.
        filenames = [
            local.path(directory) / (name + suffix)
            for directory, name, suffix in zip(dirs, names, suffixes)
        ]
        output_dict = dict(zip(ids, filenames))
        if properties:
            output_dict.update(properties)

        # Generate a partial function with meta keyword arg predefined.
        partial_wrapped_func = partial(wrapped, meta=output_dict)

        def replacement_fixture_func(**kwargs):
            # Function which will replace the decorated function. Run the
            # wrapped function with injected `meta` variable. Then return the
            # metadata.
            partial_wrapped_func(**kwargs)
            return output_dict

        wrapped_argnames = get_param_names(wrapped)

        # pytest checks that all fixture arguments are valid fixtures, so we
        # need to purge all references to `meta` parameter
        wrapped_argnames.remove("meta")
        del kwargs["meta"]
        replacement_fixture_func = replace_param_sig(
            replacement_fixture_func, wrapped_argnames
        )
        return replacement_fixture_func(**kwargs)

    return wrapper(wrapped)
