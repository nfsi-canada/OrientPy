import os.path
from os import listdir
import re
from numpy.distutils.core import setup
from pathlib import Path

def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname) as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")

scripts = [str(x) for x in Path('Scripts').iterdir() if x.is_dir()]

setup(
    name='orientpy',
    version=find_version('orientpy', '__init__.py'),
    description='Python Module for computing station orientations',
    author='Pascal Audet',
    maintainer='Pascal Audet',
    maintainer_email='pascal.audet@uottawa.ca',
    classifiers=[
         'Development Status :: 3 - Alpha',
         'License :: OSI Approved :: MIT License',
         'Programming Language :: Python :: 3.6',
         'Programming Language :: Python :: 3.7'],
    install_requires = ['numpy', 'obspy', 'stdb', 'geographiclib'],
    python_requires =  '>=3.6',
    packages= ['orientpy'],
    scripts=scripts,
    url='')
