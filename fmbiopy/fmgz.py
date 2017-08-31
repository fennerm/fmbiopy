import os
from pathlib import Path

def gunzip(p, out_dir):
    """ Gunzip a file and output to out_dir. Keep the original file"""
    unzip_name = str(p.name)[:-3]
    unzip_path = Path(out_dir, unzip_name)
    os.system('gunzip -c '+ str(p) + ' > ' + str(unzip_path))
    return unzip_path

def gzip(f):
    """ Gzip a file """
    os.system('gzip ' + str(f))
    renamed = Path(str(f) + '.gz')
    return renamed
