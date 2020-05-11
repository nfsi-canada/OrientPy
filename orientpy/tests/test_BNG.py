import stdb
import numpy as np
from pkg_resources import resource_filename
from orientpy import BNG, utils
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from . import get_meta


def test_init_BNG():

    sta = get_meta.get_stdb()
    bng = BNG(sta)
    assert isinstance(bng, BNG), 'Failed initializing BNG object'
    return bng


def test_add_cat():

    bng = test_init_BNG()
    cat = get_meta.get_cat()
    for ev in [cat[0]]:
        accept = bng.add_event(
            ev, gacmin=30., gacmax=90.,
            depmax=40., returned=True)
        assert accept, 'Event not accepted'
    return bng


def test_add_data():

    bng = test_add_cat()

    # Get travel time info
    tpmodel = TauPyModel(model='iasp91')

    # Get Travel times
    arrivals = tpmodel.get_travel_times(
        distance_in_degree=bng.meta.gac,
        source_depth_in_km=bng.meta.dep,
        phase_list=['P', 'PP'])

    # Get first P wave arrival among P and PP
    arrival = arrivals[0]

    # Attributes from parameters
    bng.meta.ttime = arrival.time
    bng.meta.phase = arrival.name

    # Get data
    t1 = arrival.time - 15.
    t2 = arrival.time + 15.
    has_data = bng.download_data(
        client=Client(), stdata=[],
        ndval=0., new_sr=2., t1=t1, t2=t2, 
        returned=True, verbose=False)

    assert has_data, 'No data'
    return bng


def test_calc():

    bng = test_add_data()
    bng.calc(1., 15., [-2.,5.], showplot=False)

    assert bng.meta.phi is not None
    assert bng.meta.snr is not None
    assert bng.meta.cc is not None
    assert bng.meta.TR is not None
    assert bng.meta.RZ is not None
    return bng.meta

def test_average():

    phi = []; cc = []; snr = []; TR = []; RZ = []; baz = []; mag = []
    meta = test_calc()
    phi = np.array([meta.phi])
    snr = np.array([meta.snr])
    cc = np.array([meta.cc])
    TR = np.array([meta.TR])
    RZ = np.array([meta.RZ])
    baz = np.array([meta.baz])
    mag = np.array([meta.mag])

    # Set conditions for good result
    snrp = snr>-10.
    ccp = cc>-1.
    TRp = TR>-1.e3
    RZp = RZ>-1.e3

    # Indices where conditions are met
    ind = snrp*ccp*TRp*RZp

    # Get estimate and uncertainty
    val, err = utils.estimate(phi, ind)



