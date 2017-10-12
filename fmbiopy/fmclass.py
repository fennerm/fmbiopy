"""Function and classes involved in manipulating and inspecting objects"""

import importlib
import inspect
from typing import List
from typing import Type


def classname(cls: Type[object]) -> str:
    """Return the name of a class"""
    return cls.__name__


def list_classes(
        module_name: str,
        exclude: List[str] = None,
        of_type: List[str] = None) -> List[str]:
    """List all classes in a python module

    Parameters
    ----------
    module_name
        Module name to inspect
    exclude, Optional
        List of classes to exclude from the return
    type, Optional
        Only classes of given type should be returned

    Returns
    -------
    A list of classes in `module_name`

    """
    module = importlib.import_module(module_name)
    classes = []

    # List all classes
    for _, obj in inspect.getmembers(module):
        if inspect.isclass(obj):
            if obj.__module__ == module.__name__:
                if of_type:
                    # List object's class and class inheritance
                    class_inheritance = inspect.getmro(obj)
                    # Get class inheritance names
                    inheritance_names = [classname(cls)
                                         for cls in class_inheritance]
                    if of_type in inheritance_names:
                        classes.append(obj)
                else:
                    classes.append(obj)

    # Exclude some
    if exclude:
        classes = [cls for cls in classes if cls.__name__ not in exclude]

    return sorted(classes, key=classname)  # type: ignore
