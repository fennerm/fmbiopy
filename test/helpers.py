"""Test helper functions."""
from random import randint
from uuid import uuid4

from pandas import DataFrame
from numpy.random import binomial
from plumbum import local
from plumbum.cmd import picard

from fmbiopy.paths import is_empty


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
