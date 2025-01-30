import setuptools
import os.path
from os import listdir
import re
from setuptools import setup
from pathlib import Path


def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname) as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name='orientpy',
    version=find_version('orientpy', '__init__.py'),
    description='Python Module for computing station orientations',
    author='Pascal Audet',
    author_email='pascal.audet@uottawa.ca',
    maintainer='Pascal Audet',
    maintainer_email='pascal.audet@uottawa.ca',
    url='https://github.com/nfsi-canada/OrientPy',
    download_url='https://github.com/nfsi-canada/OrientPy/archive/refs/tags/v0.0.2.tar.gz',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'],
    install_requires=['numpy', 'obspy', 'stdb', 'geographiclib'],
    python_requires='>=3.6',
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={'orientpy': ['dispmaps/R*txt']},
    entry_points={'console_scripts':
        ['bng_calc_auto=orientpy.scripts.bng_calc_auto:main',
         'bng_average=orientpy.scripts.bng_average:main',
         'dl_calc=orientpy.scripts.dl_calc:main',
         'dl_average=orientpy.scripts.dl_average:main']})
