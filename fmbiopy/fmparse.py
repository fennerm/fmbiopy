"""Parsing utilities"""
from collections import defaultdict
import sys
from typing import (
        Any,
        Dict,
        List,
        Sequence,
        Tuple,
        Union,
        )

from docopt import (
        docopt,
        DocoptExit,
        )

def helpful_docopt(doc: str, *args, **kwargs)-> Dict[str, Any]:
    """Wrapper around docopt which prints the full usage message upon error"""
    try:
        opt = docopt(doc, help=False, *args, **kwargs)
        return opt
    except DocoptExit:
        print(doc)
        sys.exit()
