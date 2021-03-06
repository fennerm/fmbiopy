"""Test suite for fmbiopy.fmlist"""
from pytest import mark, raises

from fmbiopy.iter import *


def test_any_endswith():
    assert any_endswith(["bar", "foo"], "o")
    assert not any_endswith(["bar", "foo"], "x")


def test_all_equal():
    assert not all_equal([1, 3])
    assert all_equal(["foo", "foo"])


def test_as_strs():
    # Trivial
    as_strs([1, 2, 3])


def test_distribute_across():
    with raises(IndexError):
        distribute_across(4, [])

    assert distribute_across(2, [1, 1]) == [2, 2]
    assert distribute_across(5, [4]) == [9]
    assert distribute_across(2, [0, 0, 2]) == [1, 1, 2]


@mark.parametrize(
    "inp,expected",
    [(3, [3]), (None, [None]), ([1, 2], [1, 2]), ("foo", ["foo"])],
)
def test_ensure_list(inp, expected):
    assert ensure_list(inp) == expected


@mark.parametrize(
    "inp,expected",
    [
        (["a", "", "b"], ["a", "b"]),
        (["a", None, "b"], ["a", "b"]),
        ([1, 2], [1, 2]),
        ([1, 2, []], [1, 2]),
    ],
)
def test_exclude_blank(inp, expected):
    assert exclude_blank(inp) == expected


@mark.parametrize(
    "inp,expected",
    [
        ([[1, 2], 1], [1, 2, 1]),
        ([1, 2, 1], [1, 2, 1]),
        ([[[1], 2], 3], [1, 2, 3]),
    ],
)
def test_flatten(inp, expected):
    assert flatten(inp) == expected


@mark.parametrize(
    "inp, expected",
    [
        ([1, 1, 2], [1, 2]),
        ([2, None], [2, None]),
        ("abc", ["a", "b", "c"]),
        ([], []),
    ],
)
def test_get_unique(inp, expected):
    assert set(get_unique(inp)) == set(expected)


@mark.parametrize(
    "inp1, inp2, expected",
    [([1, 2], [3, 4], [1, 3, 2, 4]), ([1, None], [3, 4], [1, 3, None, 4])],
)
def test_interleave(inp1, inp2, expected):
    assert interleave(inp1, inp2) == expected


def test_interleave_raises_err():
    with raises(ValueError):
        interleave([1, 2], [1])


@mark.parametrize(
    "inp, expect",
    [([True, False], False), ([False], True), ([False, False], True)],
)
def test_none(inp, expect):
    assert none(inp) == expect


def test_split_into_chunks():
    assert split_into_chunks("", 1) == [""]
    assert split_into_chunks("ab", 2) == ["a", "b"]
    assert split_into_chunks("abc", 2) == ["ac", "b"]

    xs = list(range(8))
    result = split_into_chunks(xs, 3)
    assert [len(x) <= 3 for x in result]
    assert set(flatten(result)) == set(xs)


@mark.parametrize(
    "inp,expect",
    [
        ([1, 2, 0, 3, 0, 4], [[1, 2], [3], [4]]),
        ([1, 2, 0], [[1, 2]]),
        ([], [[]]),
        ([1, 2], [[1, 2]]),
    ],
)
def test_split_list(inp, expect):
    assert split_list(inp, 0) == expect


def test_split_list_err():
    with raises(TypeError):
        split_list("ab", 0)
