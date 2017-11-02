"""Custom exceptions"""

class EmptyListError(ValueError):
    """Raised when an empty list is passed invalidly"""
    pass

class ParseError(ValueError):
    """Raised when an exception occurs during file parsing"""
    pass
