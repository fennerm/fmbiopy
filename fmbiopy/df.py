"""Utilities related to pandas data frames."""


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
