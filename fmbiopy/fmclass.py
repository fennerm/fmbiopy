"""Function and classes involved in manipulating and inspecting objects"""

import importlib
import inspect
import typing


def list_classes(module_name: str, exclude=None) -> typing.List:
    """List all classes in a python module

    Parameters
    ----------
    module_name
        Module name to inspect
    exclude : List[class], Optional
        List of classes to exclude from the return

    Returns
    -------
    A list of classes in `module_name`

    """
    module = importlib.import_module(module_name)
    tasks = []

    # List all classes
    for name, obj in inspect.getmembers(module):
        if inspect.isclass(obj):
            tasks.append(obj)

    # Exclude some
    if exclude:
        tasks = [task for task in tasks if task not in exclude]

    return tasks
