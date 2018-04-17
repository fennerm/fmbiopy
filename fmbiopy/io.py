"""Utilities related to file import/output."""
import csv

from plumbum import local


def list_to_csv(l, output_file, delimiter=','):
    """Write a list of rows to a csv or tsv."""
    output_file = local.path(output_file)
    if not output_file.dirname.exists():
        output_file.dirname.mkdir()

    with output_file.open('w') as f:
        writer = csv.writer(f, delimiter=delimiter)
        for row in l:
            writer.writerow(row)


def write_table(df, outfile, delimiter=','):
    """Write a data frame to a file with sensible defaults.

    If the output file directory doesn't exist it is created. The row index is
    not included in the output file.
    """
    if not outfile.dirname.exists():
        outfile.dirname.mkdir()
    df.to_csv(str(outfile), sep=delimiter, index=False)
