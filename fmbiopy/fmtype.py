""" Custom types for static typing with the typing module and mypy"""

from pathlib import (
        Path,
        PurePath,
        )
from typing import (
        TypeVar,
        Union,
        )

"""Generic TypeVar for specifying that output type is based upon input type"""
T = TypeVar('T')

StrOrBytes = Union[str, bytes]

PathLike = Union[str, Path, PurePath]
