"""Test suite for fmbiopy.df."""
from pandas import DataFrame
from pandas.errors import MergeError
from pytest import fixture, lazy_fixture, mark, raises

from fmbiopy.df import *
from test.helpers import assert_df_equals


@fixture
def grouped_df():
    data = {"a": [1, 1, 1], "b": [2, 2, 2], "group": ["g1", "g1", "g2"]}
    return DataFrame.from_dict(data)


@fixture
def eg_df():
    return DataFrame.from_dict(
        {"a": ["foo", "bar", "car"], "b": ["oof", "rab", "rac"], "c": [1, 2, 3]}
    )


@fixture
def empty_df():
    return DataFrame(columns=["a", "b", "c"])


@mark.parametrize("include_by", [False, True])
def test_split(grouped_df, include_by):
    dfs = split(grouped_df, by="group", include_by=include_by)
    if include_by:
        expected_ncol = 3
    else:
        expected_ncol = 2

    assert len(dfs) == 2
    assert all([len(df.columns) == expected_ncol for df in dfs])
    assert dfs[0].shape[0] == 2
    assert dfs[1].shape[0] == 1


@mark.parametrize(
    "df1,df2,output",
    [
        (
            lazy_fixture("eg_df"),
            DataFrame({"a": ["foo", "dad", "car"], "b": ["oof", "dad", "rar"]}),
            {"a": ["bar", "car"], "b": ["rab", "rac"], "c": [2, 3]},
        ),
        (
            lazy_fixture("empty_df"),
            lazy_fixture("eg_df"),
            lazy_fixture("empty_df"),
        ),
        (
            lazy_fixture("eg_df"),
            lazy_fixture("empty_df"),
            lazy_fixture("eg_df"),
        ),
    ],
)
def test_df_subtract(df1, df2, output):
    assert_df_equals(df_subtract(df1, df2), output)


def test_df_subtract_no_matching_cols_raises_err(eg_df):
    df2 = eg_df.rename(columns={"a": "x", "b": "y", "c": "z"})
    with raises(MergeError):
        df_subtract(eg_df, df2)
