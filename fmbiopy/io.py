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


def write_table(df, outfile, **kwargs):
    """Write a data frame to a file with sensible defaults.

    If the output file directory doesn't exist it is created. The row index is
    not included in the output file.
    """
    outfile = local.path(outfile)
    if not outfile.dirname.exists():
        outfile.dirname.mkdir()
    df.to_csv(str(outfile), index=False, **kwargs)


def write_table_with_header(df, header, filename, **kwargs):
    """Write a table with a prepended header."""
    filename = local.path(filename)
    if not filename.dirname.exists():
        filename.dirname.mkdir()
    with filename.open('a') as f:
        f.writelines(header)
        df.to_csv(f, index=False)


def read_header(filename, comment_char='#'):
    """Read the header from a file

    Returns
    -------
    List[str]

    """
    header = []
    with open(filename, 'r') as f:
        while True:
            line = f.readline()
            if line.startswith(comment_char) or line.isspace():
                header.append(line)
            else:
                break
    return header
