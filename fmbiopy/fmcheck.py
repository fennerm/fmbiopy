import os

""" Functions which start with check_ raise exceptions, all others return
Bools """

def check_non_empty(l):
    """ Raise exception if list is empty """

    if not l:
        raise ValueError("List is empty")
    else:
        return True

def all_equal(l):
    """ Test whether all items in list are equal """

    check_non_empty(l)
    return all(x == l[0] for x in l)

def any_dont_exist(paths):
    """ Return True if any path in list does not exist """

    check_non_empty(paths)
    exists = map(os.path.exists, paths)
    return any(not x for x in exists)

def check_all_exist(paths):
    """ Raise an exception if any paths in list do not exist """

    if any_dont_exist(paths):
        raise ValueError("Path doesn't exist: ")
    else:
        return True

def check_suffix(x, suffixes):                                                  
    """ Check if a string x ends with any of a list of suffixes """

    if isinstance(suffixes, basestring):
        suffixes = [suffixes]

    correct = any(x.endswith(suffix) for suffix in suffixes)  
    if not correct:
        raise ValueError(x + " does not have the correct suffix " + 
                         ' '.join(suffixes))
    else:
        return True

def check_all_suffix(xs, suffixes):                                             
    """ Check if all strings in a list end with one of a list of suffixes. 
        If not, raise an exception """

    if isinstance(suffixes, basestring):
        suffixes = [suffixes]

    all_correct = all([check_suffix(x, suffixes) for x in xs])
    if not all_correct:
        raise ValueError("Incorrect suffix " + ' '.join(xs))
    else:
        return True
  

