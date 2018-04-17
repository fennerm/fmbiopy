"""Test suite for fmbiopy.df."""
import pandas as pd
from pytest import (
    fixture,
    mark,
)

from fmbiopy.df import *


@fixture
def grouped_df():
    data = {
        'a': [1, 1, 1],
        'b': [2, 2, 2],
        'group': ['g1', 'g1', 'g2']
    }
    return pd.DataFrame.from_dict(data)


@mark.parametrize('include_by', [False, True])
def test_split(grouped_df, include_by):
    dfs = split(grouped_df, by='group', include_by=include_by)
    if include_by:
        expected_ncol = 3
    else:
        expected_ncol = 2

    assert len(dfs) == 2
    assert all([len(df.columns) == expected_ncol for df in dfs])
    assert dfs[0].shape[0] == 2
    assert dfs[1].shape[0] == 1
