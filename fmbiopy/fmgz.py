import os

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

