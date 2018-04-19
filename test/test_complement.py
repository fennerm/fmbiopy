"""Test suite for exclude_shared_rows.py."""
from uuid import uuid4

import pandas as pd
from plumbum import local
from pytest import (
    fixture,
    mark,
)

from bin.complement import complement
from fmbiopy.io import list_to_csv


def gen_data():
    """
    Generate a list of table rows in which the 1st, 2nd rows are always
    identical, the 3rd row is always unique and the 4th row is identical to the
    1st except in the 3rd column.
    """
    data = [
        ['n1', 'n2', 'n3'],
        ['1', 'a', 'b'],
        ['2', 'c', 'd'],
        [uuid4().hex for _ in range(3)],
        ['1', 'a', uuid4().hex]
    ]
    return data


@fixture(params=[',', '\t'])
def tables(request, sandbox):
    delimiter = request.param
    nsamples = 3
    if delimiter == ',':
        extension = '.csv'
    else:
        extension = '.tsv'
    table_files = [sandbox / (uuid4().hex + extension) for _ in range(nsamples)]
    datatables = [gen_data() for _ in range(nsamples)]
    for table, file in zip(datatables, table_files):
        list_to_csv(table, file, delimiter=delimiter)

    output_dict = {
        'data': datatables,
        'files': table_files,
        'delimiter': delimiter,
    }
    yield output_dict
    for f in table_files:
        f.delete()


@mark.parametrize('include_cols', [['0', '1'], None])
def test_complement(sandbox, tables, include_cols):
    output_prefix = sandbox + '/' + uuid4().hex + '/'
    complement(tables['files'],
               output_prefix=output_prefix,
               delimiter=tables['delimiter'],
               include_cols=include_cols)

    output_files = [local.path(output_prefix) / f.name for f in tables['files']]
    output_dfs = [pd.read_csv(f, sep=tables['delimiter']) for f in output_files]

    if include_cols:
        expected_nrows = 1
    else:
        expected_nrows = 2

    assert all([df.shape[0] == expected_nrows for df in output_dfs])


@mark.parametrize('output_prefix', [uuid4().hex, uuid4().hex + '/'])
def test_output_prefix_behavior(sandbox, output_prefix, tables):
    output_prefix = sandbox + '/' + output_prefix
    expected_output = [local.path(output_prefix + f.name)
                       for f in tables['files']]
    complement(tables['files'],
               output_prefix=output_prefix,
               delimiter=tables['delimiter'],
               include_cols=None)
    for f in expected_output:
        assert f.exists()
    for f in expected_output:
        f.delete()


def test_target_param(sandbox, tables):
    output_prefix = sandbox + '/' + uuid4().hex + '/'
    target = local.path(output_prefix + tables['files'][0].name)
    nontargets = [local.path(output_prefix + f.name)
                  for f in tables['files'][1:]]
    complement(tables['files'],
               output_prefix=output_prefix,
               delimiter=tables['delimiter'],
               include_cols=None,
               target=target)
    assert target.exists()
    for f in nontargets:
        assert not f.exists()
