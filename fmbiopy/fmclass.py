"""Function and classes involved in manipulating and inspecting objects"""

from importlib import import_module
from inspect import (
        getmembers,
        getmro,
        isclass,
        )

from typing import (
        Any,
        List,
        Tuple,
        Type,
        )


def classname(cls: Type) -> str:
    """Return the name of a class"""
    return cls.__name__


def list_classes(
        module: str,
        package: str = None,
        exclude: List[str] = None,
        of_type: List[str] = None) -> List[Any]:
    """List all classes in a python module

    Parameters
    ----------
    module_name
        Module name to inspect
    package: optional
        If doing a relative import, specify the package name to import from
    exclude, optional
        List of classes to exclude from the return
    type, optional
        Only classes of given type should be returned

    Returns
    -------
    A list of classes in `module_name`

    """
    import_module(package)
    imported = import_module(''.join(['.', module]), package)
    classes: List[Type] = []

    # List all classes
    for _, obj in getmembers(imported):
        if isclass(obj):
            if of_type:
                if obj.__module__ == imported.__name__:
                    # List object's class and class inheritance
                    class_inheritance: Tuple[Type, ...] = getmro(obj)

                    # Get class inheritance names
                    inheritance_names: List[str] = [
                            classname(cls) for cls in class_inheritance]

                    for typ in of_type:
                        if typ in inheritance_names:
                            classes.append(obj)
            else:
                classes.append(obj)

    # Exclude some
    if exclude:
        classes = [cls for cls in classes if classname(cls) not in exclude]

    return sorted(classes, key=classname)

def parent_names(cls: Type):
    """Get the names of the parent classes of `cls`"""
    parent_classes = cls.__bases__
    names = [classname(parent) for parent in parent_classes]
    return names
