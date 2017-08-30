from setuptools import setup, find_packages

setup(name='fmbiopy',
      version='0.1',
      description="""Misc. personal bioinformatics python modules which I share 
                     between projects""",
      url='http://github.com/fennerm/fmbiopy',
      author='Fenner Macrae',
      author_email='fennermacrae@gmail.com',
      license='MIT',
      packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", 
                                      "tests"]),
      zip_safe=False)
