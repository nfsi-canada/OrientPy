import stdb
import numpy as np
import shutil
from pkg_resources import resource_filename
from pathlib import Path


dbfile = resource_filename('orientpy',
                           'examples/data/LOBS3.pkl')

curdir = Path.cwd()
bngdir = curdir / 'BNG_RESULTS'

def test_01_bng_calc():
    from orientpy.scripts import bng_calc_auto as bng
    args0 = bng.get_bng_calc_arguments(
        [dbfile, '--times', ' -5.,15.',
        '--start', '2014-10-01', '--end', '2014-12-01',
        '--window', '60.', '--bp', '0.04,0.1',
        '--min-mag', '6.', '--min-dist', '30.'])
    print('*'*20)
    print(args0)
    print('*'*20)
    bng.main(args=args0)
def test_02_bng_average():
    from orientpy.scripts import bng_average as bng
    args0 = bng.get_bng_average_arguments(
        [dbfile, '--keys', 'LOBS3', '--plot', '--save'])
    bng.main(args=args0)

def test_03_rmtree():
    shutil.rmtree(bngdir)
