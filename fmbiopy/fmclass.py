"""Function and classes involved in manipulating and inspecting objects"""

import importlib
import inspect
from typing import List
from typing import Type


def classname(cls: Type[object]) -> str:
    """Return the name of a class"""
    return cls.__name__


def list_classes(
        module: str,
        package: str = None,
        exclude: List[str] = None,
        of_type: List[str] = None) -> List[str]:
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
    importlib.import_module(package)
    imported = importlib.import_module(''.join(['.', module]), package)
    classes = []

    # List all classes
    for _, obj in inspect.getmembers(imported):
        if inspect.isclass(obj):
            if obj.__module__ == imported.__name__:
                if of_type:
                    # List object's class and class inheritance
                    class_inheritance = inspect.getmro(obj)
                    # Get class inheritance names
                    inheritance_names = [classname(cls)
                                         for cls in class_inheritance]
                    for typ in of_type:
                        if typ in inheritance_names:
                            classes.append(obj)
                else:
                    classes.append(obj)

    # Exclude some
    if exclude:
        classes = [cls for cls in classes if cls.__name__ not in exclude]

    return sorted(classes, key=classname)  # type: ignore
