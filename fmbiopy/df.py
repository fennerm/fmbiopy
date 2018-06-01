"""Utilities related to pandas data frames."""
from pandas import merge


def get_colnames(df):
    """Get the column names of a data frame as a list."""
    # Faster than using list(df)
    return df.columns.get_values().tolist()


def split(df, by, include_by=True):
    """Split a data frame by a grouping column.

    Parameters
    ----------
    df: pandas.DataFrame
        The input data frame.
    by: str
        The name of the column to split by.
    include_by: bool
        If True, the group column is conserved in the output tables, if False it
        is removed.

    Returns
    -------
    List[pandas.DataFrame]

    """
    split_dfs = [rows for _, rows in df.groupby(by)]
    if not include_by:
        split_dfs = [x.drop(by, axis=1).reindex() for x in split_dfs]
    return split_dfs


def df_subtract(df1, df2):
    """Remove rows in `df` which are also in `df2` (both DataFrames)."""
    merged = merge(df1, df2, how="outer", indicator=True)
    return (
        merged[merged["_merge"] == "left_only"]
        .drop("_merge", axis=1)
        .reset_index(drop=True)
    )
