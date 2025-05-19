
import numpy as np
from obspy.clients.fdsn import Client
from obspy import UTCDateTime, read, Stream, Trace
import os
import math
import copy
from numpy import nan, isnan, abs


def floor_decimal(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n * multiplier) / multiplier


def traceshift(trace, tt):
    """
    Function to shift traces in time given travel time

    Parameters
    ----------
    trace : :class:`~obspy.core.Trace`
        Seismogram 
    tt : float
        Time delay (sec)

    Returns
    -------
    rtrace : :class:`~obspy.core.Trace`
        Shifted copy of the input trace

    """

    # Define frequencies
    nt = trace.stats.npts
    dt = trace.stats.delta
    freq = np.fft.fftfreq(nt, d=dt)

    # Fourier transform
    ftrace = np.fft.fft(trace.data)

    # Shift
    for i in range(len(freq)):
        ftrace[i] = ftrace[i]*np.exp(-2.*np.pi*1j*freq[i]*tt)

    # Back Fourier transform and return as trace
    rtrace = trace.copy()
    rtrace.data = np.real(np.fft.ifft(ftrace))

    # Update start time
    rtrace.stats.starttime -= tt

    return rtrace


def download_data(client=None, sta=None, start=UTCDateTime(),
                  end=UTCDateTime(), verbose=False,
                  zcomp='Z'):
    """
    Function to build a stream object for a seismogram in a given time window either
    by downloading data from the client object or alternatively first checking if the
    given data is already available locally.

    Parameters
    ----------
    client : :class:`~obspy.client.fdsn.Client`
        Client object
    sta : :class:`~stdb.StDbElement`
        Station metadata from :mod:`~StDb` data base
    start : :class:`~obspy.core.utcdatetime.UTCDateTime`
        Start time for request
    end : :class:`~obspy.core.utcdatetime.UTCDateTime`
        End time for request
    verbose : bool
        Whether or not to print messages to screen during run-time
    zcomp : str
        Vertical Component Identifier. Should be a single character.
        This is different then 'Z' only for fully unknown component
        orientation (i.e., components are 1, 2, 3)

    Returns
    -------
    err : bool
        Boolean for error handling (`False` is associated with success)
    trN : :class:`~obspy.core.Trace`
        Trace of North component of motion
    trE : :class:`~obspy.core.Trace`
        Trace of East component of motion
    trZ : :class:`~obspy.core.Trace`
        Trace of Vertical component of motion

    """

    from fnmatch import filter
    from obspy import read, Stream
    from os.path import dirname, join, exists
    from numpy import any
    from math import floor

    for loc in sta.location:

        # Construct location name
        if loc == "--":
            tloc = ""
        else:
            tloc = copy.copy(loc)

        # Construct Channel List
        cha = sta.channel.upper() + '?'
        msg = "*          {0:s}.{1:2s}?.{2:2s} - Checking Network".format(
            sta.station, sta.channel.upper(), tloc)
        print(msg)

        # Get waveforms, with extra 1 second to avoid
        # traces cropped too short - traces are trimmed later
        try:
            st = client.get_waveforms(
                network=sta.network,
                station=sta.station,
                location=tloc,
                channel=cha,
                starttime=start,
                endtime=end+1.)
        except Exception:
            print("*              - No Data Available")
            st = None

        # check if download got all needed data
        if st is not None and len(st.select(component=zcomp)) >= 1 and len(st.select(component="N")) >= 1 and len(st.select(component='E')) >= 1:
            print("*              - " + zcomp.upper() + "NE Data Downloaded")
            break

        # check if download got all needed data
        elif st is not None and len(st.select(component=zcomp)) >= 1 and len(st.select(component="1")) >= 1 and len(st.select(component='2')) >= 1:
            print("*              - " + zcomp.upper() + "12 Data Downloaded")
            break
        else:
            print("*              - Stream is missing components")

    # Check the correct 3 components exist
    if st is None:
        print("* Error retrieving waveforms")
        print("**************************************************")
        return True, None

    # Three components successfully retrieved
    else:

        st.merge()

        # Don't want any masked data or data with nans
        data_error = [np.ma.count_masked(tr.data) > 0 for tr in st]
        if np.sum(data_error) > 0:
            print("* Missing Data Error")
            print("* -> Aborting")
            print("**************************************************")
            return True, None

        # Detrend and apply taper
        st.detrend().detrend('linear').taper(
            max_percentage=0.05, max_length=5.)

        # Check start times
        if not np.all([tr.stats.starttime == start for tr in st]):
            print("* Start times are not all close to true start: ")
            [print("*   "+tr.stats.channel+" " +
                   str(tr.stats.starttime)+" " +
                   str(tr.stats.endtime)) for tr in st]
            print("*   True start: "+str(start))
            print("*        -> Shifting traces to true start")
            delay = [tr.stats.starttime - start for tr in st]
            st_shifted = Stream(
                traces=[traceshift(tr, dt) for tr, dt in zip(st, delay)])
            st = st_shifted.copy()

        # Check sampling rate
        sr = st[0].stats.sampling_rate
        sr_round = float(floor_decimal(sr, 0))
        if not sr == sr_round:
            print("* Sampling rate is not an integer value: ", sr)
            print("* -> Resampling")
            st.resample(sr_round, no_filter=False)

        # Try trimming
        try:
            st.trim(start, end)
        except Exception:
            print("* Unable to trim")
            print("* -> Aborting")
            print("**************************************************")
            return True, None

        # Check final lengths - they should all be equal if start times 
        # and sampling rates are all equal and traces have been trimmed
        if not np.allclose([tr.stats.npts for tr in st[1:]], st[0].stats.npts):
            print("* Lengths are incompatible: ")
            [print("*     "+str(tr.stats.npts)) for tr in st]
            print("*     Trimming to shortest segment")
            L = int(np.unique(np.min(np.array([len(tr.data) for tr in st]))))
            st[0].data = st[0].data[0:L]
            st[1].data = st[1].data[0:L]
            st[2].data = st[2].data[0:L]

        print("* Waveforms Retrieved...")
        return False, st
