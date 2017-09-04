import os
from pathlib import Path

def gunzip(path: Path, out_dir: str) -> Path:
    """ Gunzip a file and output to out_dir. Keep the original file"""
    unzip_name = str(path.name)[:-3]
    unzip_path = Path(out_dir, unzip_name)
    os.system('gunzip -c '+ str(path) + ' > ' + str(unzip_path))
    return unzip_path

def gzip(path: Path) -> Path:
    """ Gzip a file """
    os.system('gzip ' + str(path))
    renamed = Path(str(path) + '.gz')
    return renamed
