"""Utilities for working with readcount files.

At present most functions are designed to work with output of IARCbioinfo's
fork of mpileup2readcounts (github.com/IARCbioinfo/mpileup2readcounts), but
other data formats might be added later.

IMPORTANT: loc column in mpileup2readcounts output must be renamed to pos so as
not to conflict with loc based indexing.
"""
from collections import OrderedDict

from numpy import NaN
from pandas import concat, isnull

from fmbiopy.fmlist import pairwise_intersect


def combine_indel_counts(indel_cell):
    """Combine indel counts from a cell in a count table."""
    if isnull(indel_cell):
        return NaN
    indel_counts = OrderedDict()
    for indel in indel_cell.split("|"):
        indel_count, indel_type = tuple(indel.split(":"))
        try:
            indel_counts[indel_type.upper()] += int(indel_count)
        except KeyError:
            indel_counts[indel_type.upper()] = int(indel_count)
    return "|".join([":".join([str(v), k]) for k, v in indel_counts.items()])


def destrand_counts(counts):
    """Combine counts in table across strands.

    Parameters
    ----------
    counts: DataFrame
        mpileup2readcounts output as a pandas DataFrame

    Returns
    -------
    DataFrame

    """
    combined_counts = counts[["chr", "pos", "ref", "depth"]].copy()
    for i, nuc in enumerate(["A", "T", "C", "G"]):
        combined_counts[nuc] = counts[nuc] + counts[nuc.lower()]
    for indel_type in ["Insertion", "Deletion"]:
        combined_counts[indel_type] = counts[indel_type].apply(
            combine_indel_counts
        )
    return combined_counts


def find_shared_bases(per_sample_counts, nucleotide):
    nsamples = len(per_sample_counts)
    base_counts = concat(
        [counts[nucleotide] for counts in per_sample_counts], axis=1
    )
    return ((base_counts == 0).sum(axis=1) < (nsamples - 1)).values


def _get_indel_set(indel_cell):
    """Get the indel alleles present in a readcount cell as a set."""
    if isnull(indel_cell):
        indel_set = set()
    else:
        indel_set = set(
            [indel.split(":")[1] for indel in indel_cell.split("|")]
        )
    return indel_set


def _get_shared_indels_from_row(row):
    """Get indels which are present in multiple samples from row of indel
    readcount table.
    """
    indel_sets = [_get_indel_set(cell) for cell in row.values]
    shared_indels = pairwise_intersect(indel_sets)
    if shared_indels:
        return list(shared_indels)
    else:
        return NaN


def find_shared_indels(per_sample_counts, allele):
    """Find all indel variants which are not unique to a sample.

    Parameters
    ----------
    per_sample_counts: List[DataFrame]
        List of unstranded count tables. One for each sample.
    allele: str
        One of ["Insertion", "Deletion"]

    Returns
    -------
    DataFrame
        Return a DataFrame containing shared indels and their indices.

    """
    indel_counts = concat(
        [counts[allele] for counts in per_sample_counts], axis=1
    )
    mult_samples_have_indels = indel_counts.notnull().sum(axis=1) > 1
    mult_indel_counts = indel_counts[mult_samples_have_indels]
    shared_indels = mult_indel_counts.apply(
        _get_shared_indels_from_row, axis=1, result_type="reduce"
    ).dropna()
    return shared_indels
