"""Metaprogramming functions for introspection and manipulation of objects."""
from copy import copy
from importlib import import_module
from inspect import getmembers, getmro, isclass, Parameter, signature


def get_param_names(func):
    """Get the argument names of a function."""
    return list(func.__code__.co_varnames[: func.__code__.co_argcount])


def replace_param_sig(func, params):
    """Replace the parameter signature of a function.

    Parameters
    ----------
    func: Callable
    params: Iterable[str]
        New list of argument names for the output function.

    Returns
    -------
    Callable
        A new function which is a copy of `func` with its parameter signature
        replaced.

    """
    func_copy = copy(func)
    new_parameters = [
        Parameter(param, Parameter.POSITIONAL_OR_KEYWORD) for param in params
    ]
    new_signature = signature(func_copy).replace(parameters=new_parameters)
    func_copy.__signature__ = new_signature
    return func_copy


def classname(cls):
    """Return the name of a class"""
    return cls.__name__


def list_classes(module, package=None, exclude=None, of_type=None):
    """List all classes in a python module

    Parameters
    ----------
    module: str
        Module name to inspect
    package: str, Optional
        If doing a relative import, specify the package name to import from
    exclude: List[str], optional
        List of classes to exclude from the return
    of_type: List[str], optional
        Only classes of given type should be returned

    Returns
    -------
    List[class]
        A list of classes in `module`
    """
    if package is None:
        imported = import_module(module)
    else:
        import_module(package)
        imported = import_module("".join([".", module]), package)

    classes = []

    # List all classes
    for _, obj in getmembers(imported):
        # Check if class
        if isclass(obj):
            # Check if class is defined in the target file.
            if obj.__module__ == imported.__name__:
                if of_type:
                    # List object's class and class inheritance
                    class_inheritance = getmro(obj)

                    # Get class inheritance names
                    inheritance_names = [
                        classname(cls) for cls in class_inheritance
                    ]

                    # Check if any of the objects inheritance is of the target
                    # type.
                    if any([typ in inheritance_names for typ in of_type]):
                        classes.append(obj)
                else:
                    classes.append(obj)

    # Exclude some
    if exclude:
        classes = [cls for cls in classes if classname(cls) not in exclude]

    return sorted(classes, key=classname)


def parent_names(cls):
    """Get the names of the parent classes of a class `cls`"""
    parent_classes = cls.__bases__
    names = [classname(parent) for parent in parent_classes]
    return names
