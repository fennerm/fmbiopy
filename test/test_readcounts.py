from collections import OrderedDict

import numpy as np
from pandas import DataFrame
from pytest import fixture, mark

from fmbiopy.readcounts import *


@fixture
def stranded_counts():
    return DataFrame(
        OrderedDict(
            (
                ("chr", ["scaf1", "scaf2"]),
                ("pos", [1, 2]),
                ("ref", ["G", "A"]),
                ("depth", [16, 14]),
                ("A", [1, 0]),
                ("T", [0, 2]),
                ("C", [2, 0]),
                ("G", [3, 1]),
                ("a", [1, 0]),
                ("t", [2, 3]),
                ("c", [1, 0]),
                ("g", [1, 2]),
                ("Insertion", ["2:AA|1:aa", "1:C|2:T"]),
                ("Deletion", [np.NaN, "1:AC|2:ac"]),
            )
        )
    )


@fixture
def unstranded_counts():
    return DataFrame(
        OrderedDict(
            (
                ("chr", ["scaf1", "scaf2"]),
                ("pos", [1, 2]),
                ("ref", ["G", "A"]),
                ("depth", [16, 14]),
                ("A", [2, 0]),
                ("T", [2, 5]),
                ("C", [3, 0]),
                ("G", [4, 3]),
                ("Insertion", ["3:AA", "1:C|2:T"]),
                ("Deletion", [np.NaN, "3:AC"]),
            )
        )
    )


@fixture
def empty_stranded_counts(stranded_counts):
    return DataFrame(columns=stranded_counts.columns.values)


@fixture
def empty_unstranded_counts(unstranded_counts):
    return DataFrame(columns=unstranded_counts.columns.values)


def test_destrand_counts_doesnt_raise_error_for_empty_df(
    empty_stranded_counts, empty_unstranded_counts
):
    assert destrand_counts(empty_stranded_counts).equals(
        empty_unstranded_counts
    )


def test_destrand_counts(stranded_counts, unstranded_counts):
    assert destrand_counts(stranded_counts).equals(unstranded_counts)


@mark.parametrize(
    "input,expected",
    [
        ("2:a", "2:A"),
        ("2:A|2:a", "4:A"),
        (np.NaN, np.NaN),
        ("2:AC|1:A", "2:AC|1:A"),
    ],
)
def test_combine_indel_counts(input, expected):
    np.testing.assert_equal(combine_indel_counts(input), expected)


def test_find_shared_bases_with_one_below_cutoff():
    counts = [
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [10]),
                    ("T", [10]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [18]),
                    ("T", [2]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
    ]
    assert list(find_shared_bases(counts, "T")) == [True]


def test_find_shared_bases_with_all_above_cutoff():
    counts = [
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [10]),
                    ("T", [10]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [10]),
                    ("T", [10]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
    ]
    assert list(find_shared_bases(counts, "T")) == [True]


def test_find_shared_bases_with_unique_allele():
    counts = [
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [10]),
                    ("T", [10]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [20]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
    ]
    assert list(find_shared_bases(counts, "T")) == [False]


def test_find_shared_bases_with_all_zero():
    counts = [
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [0]),
                    ("A", [0]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [0]),
                    ("A", [0]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
    ]
    assert list(find_shared_bases(counts, "T")) == [False]


def test_find_shared_indels_with_one_below_cutoff():
    counts = [
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [10]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", ["10:A"]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [19]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", ["1:A"]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
    ]
    assert list(find_shared_indels(counts, "Insertion")) == [["A"]]


def test_find_shared_indels_with_all_above_cutoff():
    counts = [
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [10]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", ["10:A"]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [10]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", ["10:A"]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
    ]
    assert find_shared_indels(counts, "Insertion").to_dict() == {0: ["A"]}


def test_find_shared_indels_with_unique_allele():
    counts = [
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [10]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", ["10:A"]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [20]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
    ]
    assert find_shared_indels(counts, "Insertion").to_dict() == dict()


def test_find_shared_indels_with_all_zero():
    counts = [
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [0]),
                    ("A", [0]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [0]),
                    ("A", [0]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
    ]
    assert find_shared_indels(counts, "Insertion").to_dict() == dict()


def test_find_shared_indels_with_multiple_shared():
    counts = [
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [0]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", ["10:A|10:AA"]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [0]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", ["10:A|10:AA"]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
    ]
    assert set(find_shared_indels(counts, "Insertion").to_dict()[0]) == set(
        ["AA", "A"]
    )


def test_find_shared_indels_with_multiple_unique():
    counts = [
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [0]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", ["10:A|10:AA"]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
        DataFrame(
            OrderedDict(
                (
                    ("chr", ["scaf1"]),
                    ("pos", [1]),
                    ("ref", ["A"]),
                    ("depth", [20]),
                    ("A", [20]),
                    ("T", [0]),
                    ("C", [0]),
                    ("G", [0]),
                    ("Insertion", [np.NaN]),
                    ("Deletion", [np.NaN]),
                )
            )
        ),
    ]
    assert find_shared_indels(counts, "Insertion").to_dict() == dict()
