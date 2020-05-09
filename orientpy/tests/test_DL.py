import stdb
import numpy as np
from pkg_resources import resource_filename
from orientpy import DL, arguments, utils
from pathlib import Path

def get_stdb():
    dbfile = resource_filename('orientpy',
                                'examples/data/LOBS3.pkl')
    db = stdb.io.load_db(dbfile)   
    return db['LOBS3']

def test_gvmap():

    gvpath = Path('../dispmaps')
    map10 = np.loadtxt(gvpath / 'R.gv.10.txt')
    map15 = np.loadtxt(gvpath / 'R.gv.15.txt')
    map20 = np.loadtxt(gvpath / 'R.gv.20.txt')
    map25 = np.loadtxt(gvpath / 'R.gv.25.txt')
    map30 = np.loadtxt(gvpath / 'R.gv.30.txt')
    map35 = np.loadtxt(gvpath / 'R.gv.35.txt')
    map40 = np.loadtxt(gvpath / 'R.gv.40.txt')
    assert map10 is not None, 'Failed! error in loading map10 disp map'
    assert map15 is not None, 'Failed! error in loading map15 disp map'
    assert map20 is not None, 'Failed! error in loading map20 disp map'
    assert map25 is not None, 'Failed! error in loading map25 disp map'
    assert map30 is not None, 'Failed! error in loading map30 disp map'
    assert map35 is not None, 'Failed! error in loading map35 disp map'
    assert map40 is not None, 'Failed! error in loading map40 disp map'

