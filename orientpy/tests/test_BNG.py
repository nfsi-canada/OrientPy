import stdb
import numpy as np
from pkg_resources import resource_filename
from orientpy import BNG, utils, plotting
from obspy.clients.fdsn import Client
from obspy.signal.rotate import rotate_rt_ne, rotate_ne_rt
from obspy.taup import TauPyModel
from obspy import Stream
from obspy.core.event.event import Event
from . import get_meta


def test_init_BNG():

    sta = get_meta.get_stdb()
    bng = BNG(sta)
    assert isinstance(bng, BNG), 'Failed initializing BNG object'
    return bng


def test_add_event(ind=0):

    bng = test_init_BNG()
    cat = get_meta.get_cat()
    for ev in [cat[ind]]:
        accept = bng.add_event(
            ev, gacmin=30., gacmax=90.,
            depmax=40., returned=True)
        assert isinstance(ev, Event)
    return bng


def test_add_data(ind=0):

    bng = test_add_event(ind)

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

    # assert has_data, 'No data'
    return has_data, bng


def test_calc(ind=0):

    has_data, bng = test_add_data(ind)
    if not has_data:
        return None
    else:
        bng.calc(1., 15., [-2.,5.])
        assert bng.meta.phi is not None
        assert bng.meta.snr is not None
        assert bng.meta.cc is not None
        assert bng.meta.TR is not None
        assert bng.meta.RZ is not None
        return bng.meta

def test_bng_waveforms(ind=0):

    has_data, bng = test_add_data(ind)
    if not has_data:
        pass

    bng.calc(1., 15., [-2.,5.])
    stream = bng.data.copy()
    stream.filter('bandpass', freqmin=0.01,
        freqmax=0.04, zerophase=True)
    trN = stream.select(component='1')[0].copy()
    trE = stream.select(component='2')[0].copy()
    azim = bng.meta.phi

    N, E = rotate_rt_ne(trN.data, trE.data, azim)
    trN.data = -1.*N
    trE.data = -1.*E

    # Update stats of streams
    trN.stats.channel = trN.stats.channel[:-1] + 'N'
    trE.stats.channel = trE.stats.channel[:-1] + 'E'

    # Store corrected traces in new stream and rotate to
    # R, T using back-azimuth
    stcorr = Stream(traces=[trN, trE])
    stcorr.rotate('NE->RT', back_azimuth=bng.meta.baz)

    # Merge original and corrected streams
    st = stream + stcorr
    plot = plotting.plot_bng_waveforms(bng, st, 15., [-2.,5.])


def test_average():

    phi = []; cc = []; snr = []; TR = []; RZ = []; baz = []; mag = []
    for ind in range(20):
        meta = test_calc(ind)
        if meta is None:
            continue
        phi.append(meta.phi)
        snr.append(meta.snr)
        cc.append(np.abs(meta.cc))
        TR.append(meta.TR)
        RZ.append(meta.RZ)
        baz.append(meta.baz)
        mag.append(meta.mag)

    phi = np.array(phi)
    cc = np.array(cc)
    snr = np.array(snr)
    TR = np.array(TR)
    RZ = np.array(RZ)
    baz = np.array(baz)
    mag = np.array(mag)

    stkey = ''

    # Set conditions for good result
    snrp = snr>-1.
    ccp = cc>-1.
    TRp = TR>-1.
    RZp = RZ>-1.

    # Indices where conditions are met
    ind = snrp*ccp*TRp*RZp

    # Get estimate and uncertainty
    val, err = utils.estimate(phi, ind)

    plot = plotting.plot_bng_conditions(stkey, snr, cc, TR, RZ, ind)


    plot = plotting.plot_bng_results(stkey, phi, snr, cc, 
        TR, RZ, baz, mag, ind, val, err)




