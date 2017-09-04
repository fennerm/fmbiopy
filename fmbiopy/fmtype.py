""" Custom types for static typing with the typing module and mypy"""

from typing import Union, Sequence
from pathlib import Path

StringOrSequence = Union[str, Sequence[str]]

PathsOrStrings = Union[Sequence[str], Sequence[Path]]

