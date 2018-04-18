#!/usr/bin/env python
"""Remove shared rows from a set of csvs or tsvs."""
import click
import pandas as pd
from plumbum import local

from fmbiopy.df import (
    get_colnames,
    split,
)
from fmbiopy.io import write_table


def read_tables(csv_files, delimiter):
    """Read a list of csv/tsv files into data frames.

    Raises
    ------
    ValueError
        If the csv files don't have the same number of columns

    """
    dfs = [pd.read_csv(f, dtype=str, sep=delimiter) for f in csv_files]

    if not all([dfs[0].shape[1] == df.shape[1] for df in dfs]):
        raise ValueError('Input files don\'t have the same number of columns')
    return dfs


def remove_shared_rows(dfs, include_cols=None):
    """Remove shared rows from a list of data frames."""
    # These will need to be added back to the data frames later.
    colnames = [list(df.columns.values) for df in dfs]

    for i, df in enumerate(dfs):
        # Replace all column names with integers so that they can be merged
        df.columns = [str(j) for j in range(df.shape[1])]
        # Add a grouping column so that the data frames can be separated after
        # merging
        df['group'] = str(i)

    merged_df = pd.concat(dfs)

    if include_cols:
        merged_df.drop_duplicates(subset=include_cols, keep=False, inplace=True)
    else:
        merged_df.drop_duplicates(
            subset=get_colnames(merged_df)[0:-1], keep=False, inplace=True)

    final_dfs = split(merged_df, by='group', include_by=False)

    for cols, df in zip(colnames, final_dfs):
        df.columns = cols

    return final_dfs


def complement(csv_files, output_prefix, delimiter=',', include_cols=None):
    """Remove all shared rows from a list of csv/tsv files."""
    dfs = read_tables(csv_files, delimiter)
    dfs = remove_shared_rows(dfs, include_cols)
    output_filenames = [local.path(output_prefix + f.name) for f in csv_files]
    for df, filename in zip(dfs, output_filenames):
        write_table(df, filename, delimiter)


@click.command()
@click.option('-o', '--output_prefix', help='Output file prefix')
@click.option('-t', '--tab_separated', is_flag=True,
              help='Table is tab separated (comma separated is default)')
@click.option('-i', '--include_cols',
              help=('List of column indices to inspect when complementing. E.g '
                    '"-i 3,4" to inspect the 4rd and 5th columns (zero '
                    'indexed)'))
@click.argument('csv_file', nargs=-1)
def cli(csv_file, output_prefix, tab_separated, include_cols):
    """Remove rows in input files which are shared between multiple files."""
    if tab_separated:
        delimiter = '\t'
    else:
        delimiter = ','

    if include_cols:
        include_cols = include_cols.split(',')

    csv_files = [local.path(f) for f in csv_file]

    complement(csv_files=csv_files,
               output_prefix=local.path(output_prefix),
               delimiter=delimiter,
               include_cols=include_cols)



if __name__ == '__main__':
    cli()
