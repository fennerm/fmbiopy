"""Test suite for fmbiopy.io."""
from uuid import uuid4

import pandas as pd
from pytest import (
    fixture,
    mark,
)

from fmbiopy.io import *


@fixture
def list_data():
    return [['n1', 'n2', 'n3'], ['a', 'b', 'c'], ['1', '2', '3']]


@mark.parametrize('delimiter', [(','), ('\t')])
def test_list_to_csv(list_data, sandbox, delimiter):
    if delimiter == ',':
        suffix = '.csv'
    else:
        suffix = '.tsv'
    output_file = sandbox / (uuid4().hex + suffix)
    list_to_csv(list_data, output_file, delimiter)
    df = pd.read_csv(output_file, sep=delimiter)
    for i, row in enumerate(df.itertuples()):
        assert tuple(row)[1:] == tuple(list_data[i+1])
