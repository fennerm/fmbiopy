import os, sys
from pathlib import Path

from contextlib import contextmanager

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

from subprocess import Popen, PIPE

## Run a bash command and log the result
## Param:
##   command  List; Bash command
##   logger_id  String; Process ID for logging instance.
##   write_log  Boolean; If True use logging module to write logs, if False 
##              don't log
## Return:
##   The exit code of the command 
def run_command(command, logger_id=None, write_log=True):

    # If command is passed as a string, convert to list
    if isinstance(command, basestring):
        # Convert
        command = command.split()

    # Remove empty list items
    command = filter(None, command)

    # Run the command
    p = Popen(command, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    
    # Try outputting to logger. Ignore if cannot find a handler
    if write_log:
        try:
            logger = logging.getLogger(logger_id)

            if stdout:
                logger.info(stdout)
            if stderr:
                logger.error(stderr)
        except:
            pass

    return (p.returncode, stdout, stderr)
