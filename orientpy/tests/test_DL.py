import stdb
import numpy as np
from pkg_resources import resource_filename
from orientpy import DL, arguments, utils
from obspy.clients.fdsn import Client


def get_stdb():
    dbfile = resource_filename('orientpy',
                                'examples/data/LOBS3.pkl')
    db = stdb.io.load_db(dbfile)   
    return db['YH.LOBS3']

def get_cat():

    sta = get_stdb()
    cat_client = Client()

    tstart = sta.startdate
    tend = sta.enddate
    minmag = 6.0
    try:
        cat = cat_client.get_events(starttime=tstart, endtime=tend,
                                    minmagnitude=minmag, maxdepth=40.)
    except:
        raise(Exception("  Fatal Error: Cannot download Catalogue"))

    return cat

def test_gvmap():

    map10 = np.loadtxt(resource_filename('orientpy',
                            'dispmaps/R.gv.10.txt'))
    map15 = np.loadtxt(resource_filename('orientpy',
                            'dispmaps/R.gv.15.txt'))
    map20 = np.loadtxt(resource_filename('orientpy',
                            'dispmaps/R.gv.20.txt'))
    map25 = np.loadtxt(resource_filename('orientpy',
                            'dispmaps/R.gv.25.txt'))
    map30 = np.loadtxt(resource_filename('orientpy',
                            'dispmaps/R.gv.30.txt'))
    map35 = np.loadtxt(resource_filename('orientpy',
                            'dispmaps/R.gv.35.txt'))
    map40 = np.loadtxt(resource_filename('orientpy',
                            'dispmaps/R.gv.40.txt'))
    assert map10 is not None, 'Failed! error in loading map10 disp map'
    assert map15 is not None, 'Failed! error in loading map15 disp map'
    assert map20 is not None, 'Failed! error in loading map20 disp map'
    assert map25 is not None, 'Failed! error in loading map25 disp map'
    assert map30 is not None, 'Failed! error in loading map30 disp map'
    assert map35 is not None, 'Failed! error in loading map35 disp map'
    assert map40 is not None, 'Failed! error in loading map40 disp map'

def test_init_DL():

    sta = get_stdb()
    dl = DL(sta)
    assert isinstance(dl, DL), 'Failed initializing DL object'
    return dl

def test_add_cat():

    dl = test_init_DL()
    cat = get_cat()
    for ev in [cat[0]]:
        # Add event to object
        accept = dl.add_event(
            ev, gacmin=0., gacmax=180.,
            depmax=40., returned=True)
        assert accept, 'event not accepted'

    return dl

def test_add_data():

    dl = test_add_cat()
    t1 = 0.
    t2 = 4.*60.*60.
    has_data = dl.download_data(
        client=Client(), stdata=[],
        ndval=0., new_sr=2., t1=t1, t2=t2, 
        returned=True, verbose=False)
    assert has_data, 'no data'
    return dl

def test_calc():

    dl = test_add_data()
    dl.calc(showplot=False)
    assert dl.meta.R1phi is not None, 'R1phi is None'
    assert dl.meta.R2phi is not None, 'R2phi is None'
    assert dl.meta.R1cc is not None, 'R1cc is None'
    assert dl.meta.R2cc is not None, 'R2cc is None'

