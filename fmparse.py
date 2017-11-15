"""Parsing utilities"""
from sys import exit

from docopt import (
        docopt,
        DocoptExit,
        )

def helpful_docopt(doc: str, *args, **kwargs):
    """Wrapper around docopt which prints the full usage message upon error"""
    try:
        docopt(doc, help=False)
    except DocoptExit:
        print(doc)
        exit()

