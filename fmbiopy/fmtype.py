""" Custom types for static typing with the typing module and mypy"""

from typing import List
from typing import Sequence
from typing import Union


StringOrSequence = Union[str, Sequence[str]]

StringOrList = Union[str, List[str]]
