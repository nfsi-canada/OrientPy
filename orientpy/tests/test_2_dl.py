import stdb
import numpy as np
import shutil
from pkg_resources import resource_filename
from pathlib import Path
import pytest


dbfile = resource_filename('orientpy',
                           'examples/data/LOBS3.pkl')

curdir = Path.cwd()
dldir = curdir / 'DL_RESULTS'

def test_01_dl_calc():
    from orientpy.scripts import dl_calc as dl
    args0 = dl.get_dl_calc_arguments(
        [dbfile, '--keys', 'LOBS3', '-O',
        '--start', '2014-10-01', '--end', '2014-10-09',
        '--verbose', '2'])
    dl.main(args=args0)
def test_02_dl_average():
    from orientpy.scripts import dl_average as dl
    args0 = dl.get_dl_average_arguments(
        [dbfile, '--keys', 'LOBS3', '--plot', '--save',
        '--cc', '0.5', '--verbose', '2'])
    dl.main(args=args0)

def test_03_rmtree():
    shutil.rmtree(dldir)
