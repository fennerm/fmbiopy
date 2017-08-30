#!/usr/bin/env python

from __future__ import print_function
from pathlib import Path
from contextlib import contextmanager
import os, sys, logging
from subprocess import Popen, PIPE

# Gunzip a given path 'p' into 'out_dir', and keep 'p'
def gunzip(p, out_dir):
	unzip_name = str(p.name)[:-3]

	if unzip_name.endswith('.fq'):
		unzip_name = unzip_name[:-2] + 'fastq'
	unzip_path = Path(out_dir, unzip_name)
	os.system('gunzip -c '+ str(p) + ' > ' + str(unzip_path))
	return unzip_path

# Gzip a file
def gzip(f):
    os.system('gzip ' + str(f))
    renamed = Path(str(f) + '.gz')
    return renamed

# Convert a list of paths to a space separated string
def paths_to_string(paths):
    s = ' '.join(str(p) for p in paths)
    return s

# Concatenate a list of files
def concat(files, out):
	as_string = paths_to_string(files)
	os.system('cat ' + as_string + ' > ' + str(out))

#Create a directory for each element of list (names) and return their paths.
def mkdirs(names, wd):
    paths = []
    for n in names:
        p = Path(wd, n)
        paths.append(p)
        mkdir(p)
    return paths

# Create directory if it doesn't already exist
def mkdir(p):
    if not Path.is_dir(p):
        Path.mkdir(p)

# Get the part of string before the first dot
def prefix(x):
    return x.split(".")[0]

# Get the part of string after the first dot
def suffix(x):
    return ''.join(x.split(".")[1:])

# Get the part of string after the last dot
def final_suffix(x):
    sp = x.split(".")
    return sp[len(sp)-1]

# Check if all items in a list are the same
def all_equal(lst):
    return all(x == lst[0] for x in lst)

# Change working directory context safely.
# Usage:
# with working_directory(dir):
# 	<code>
@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)

## Check if a string ends with any of a list of suffixes
## Param:
##   x  String
##   suffixes  List of suffixes
## Return:
##   True if string has one of the suffixes, False otherwise
def check_suffix(x, suffixes):
	return any(x.endswith(suffix) for suffix in suffixes)

## Check if a file or list of files has a valid extension. Return error if not.
## Param:
##   path  String; File path or list of file paths
##   valid_extension  String; Single valid extension or list of valid 
##                    extensions
## Return:
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

## Replace the suffix of a string
## Param:
##   x  String
##   old_suffix  Suffix to be replaced
##   new_suffix  Suffix to replace with
## Return:
##   The string with updated suffix
def replace_suffix(x, old_suffix, new_suffix):
	if not x.endswith(old_suffix):
		raise ValueError('Given suffix does not match the actual suffix (' + 
				x + ', ' + old_suffix)
	return x[:-len(old_suffix)] + new_suffix

# Set up a logging instance for a given logfile path
# https://stackoverflow.com/questions/9321741/printing-to-screen-and-writing-\
# to-a-file-at-the-same-time
def setup_log(logfile):
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s \
                                %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=logfile,
                        filemode='w')

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    # Begin logging
    logging.info("===========================================================")
    logging.info("Begin Log")
    logging.info("===========================================================")

    return console

## Run a bash command and log the result
## Param:
##   command  List; Bash command
##   logger_id  String; Process ID for logging instance.
## Return:
##   The exit code of the command 
def run_command(command, logger_id):

    p = Popen(command, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    logger = logging.getLogger(logger_id)
    
    if stdout:
        logger.info(stdout)
    if stderr:
        logger.error(stderr)

    return p.returncode

# Return True if any path in list does not exist
def any_dont_exist(paths):
    exists = map(os.path.isfile, paths)
    return any(not x for x in exists)
