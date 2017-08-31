import os

""" Functions which start with check_ raise exceptions, all others return
Bools """

def check_non_empty(l):
    """ Raise exception if list is empty """
    if not l:
        raise ValueError("List is empty")

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

## Check if a string ends with any of a list of suffixes                        
## Param-                                                                       
##   x  String                                                                  
##   suffixes  List of suffixes                                                 
## Return-                                                                      
##   True if string has one of the suffixes, False otherwise                    
def check_suffix(x, suffixes):                                                  
    return any(x.endswith(suffix) for suffix in suffixes)  

## Check if all strings in a list end with one of a list of suffixes. Raise an  
## Exception if not.                                                            
## Param                                                                        
##   xs   List of strings                                                       
##   suffixes   List of suffixes                                                
def check_all_suffix(xs, suffixes):                                             
    all_correct = all([check_suffix(x, suffixes) for x in xs])                  
    if not all_correct:                                                         
        raise ValueError("Incorrect suffix " + xs)  

## Check if a file or list of files has a valid extension. Return error if not. 
## Param-                                                                       
##   path  String; File path or list of file paths                              
##   valid_extension  Single valid extension or list of valid                   
##                    extensions                                                
## Return-                                                                      
##   True if valid extension. Error otherwise.                                  
def check_file_extension(path, valid_extension):                                
    # If given a single path convert to a single item list so that we can loop  
    # through them correctly                                                    
    if isinstance(path, basestring):                                            
        path = [path]                                                           
        valid_extension = [valid_extension]                                     
                                                                                
    for p, ext in zip(path, valid_extension):                                   
        valid = check_suffix(p, ext)                                            
        if not valid:                                                           
            raise TypeError('Invalid filetype: ' + p + ' Expected: ' + ext)     
                                                                                
    return True     
