# Copyright 2019 Pascal Audet
#
# This file is part of OrientPy.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from obspy import Trace, Stream, UTCDateTime
from scipy.stats import circmean as cmean
from scipy.stats import circstd as cstd
from scipy.stats import hmean as hm
from orientpy import io, utils, plotting
from obspy.signal.rotate import rotate_rt_ne, rotate_ne_rt
from pkg_resources import resource_filename


class Meta(object):
    """
    A Meta object contains attributes associated with the metadata
    for a single earthquake event.

    Parameters
    ----------
    time : :class:`~obspy.core.UTCDateTime`
        Origin time of earthquake
    dep : float
        Depth of hypocenter (km)
    lon : float
        Longitude coordinate of epicenter
    lat : float
        Latitude coordinate of epicenter
    mag : float
        Magnitude of earthquake
    gac : float
        Great arc circle between station and epicenter (degrees)
    epi_dist : float
        Epicentral distance between station and epicenter (km)
    baz : float
        Back-azimuth - pointing to earthquake from station (degrees)
    az : float
        Azimuth - pointing to station from earthquake (degrees)

    """

    def __init__(self, sta, event, gacmin=5., gacmax=30., depmax=1000.):

        from obspy.geodetics.base import gps2dist_azimuth as epi
        from obspy.geodetics import kilometer2degrees as k2d

        # Extract event 4D parameters
        self.time = event.origins[0].time
        self.lon = event.origins[0].longitude
        self.lat = event.origins[0].latitude
        self.dep = event.origins[0].depth

        # Check if depth is valid type
        if self.dep is not None:
            if self.dep > 1000.:
                self.dep = self.dep/1000.
        else:
            self.dep = 10.

        # Magnitude
        self.mag = event.magnitudes[0].mag
        if self.mag is None:
            self.mag = -9.

        # Calculate epicentral distance
        self.epi_dist, self.az, self.baz = epi(
            self.lat, self.lon, sta.latitude, sta.longitude)
        self.epi_dist /= 1000
        self.gac = k2d(self.epi_dist)

        if self.gac > gacmin and self.gac < gacmax and self.dep < depmax:
            self.accept = True
        else:
            self.accept = False


class Orient(object):
    """
    A Orient object contains Class attributes that associate
    station information with a single event (i.e., earthquake)
    metadata.

    Note
    ----
    The object is initialized with the ``sta`` field only, and
    other attributes are added to the object as the analysis proceeds.

    Parameters
    ----------
    sta : :class:`~stdb.StDbElement`
        Object containing station information - from :mod:`~stdb` database.
    meta : :class:`~orientpy.classes.Meta`
        Object of metadata information for single event (initially set to None)
    data : :class:`~obspy.core.Stream`
        Stream object containing the three-component seismograms
    zcomp : str
        Vertical Component Identifier. Should be a single character.
        This is different then 'Z' only for fully unknown component
        orientation (i.e., components are 1, 2, 3)

    """

    def __init__(self, sta, zcomp='Z'):

        # Attributes from parameters
        self.sta = sta

        # Initialize meta and data objects as None
        self.meta = None
        self.data = None
        self.zcomp = zcomp

    def add_event(self, event, gacmin=None, gacmax=None,
                  depmax=1000., returned=False):
        """
        Adds event metadata to Orient object.

        Parameters
        ----------
        event : :class:`~obspy.core.event`
            Event XML object
        gacmin : float
            Minimum great-arc circle distance to consider (deg)
        gacmax : float
            Maximum great-arc circle distance to consider (deg)
        depmax : float
            Maximum event depth to consider (km)
        returned : bool
            Whether or not to return the ``accept`` attribute

        Returns
        -------
        accept : bool
            Whether or not the object is accepted for further analysis

        Attributes
        ----------
        meta : :class:`~orientpy.classes.Meta`
            Meta data object

        """
        from obspy.core.event.event import Event

        if not isinstance(event, Event):
            raise Exception("Event has incorrect type")

        # Store as object attributes
        self.meta = Meta(sta=self.sta, event=event,
                         gacmin=gacmin, gacmax=gacmax,
                         depmax=depmax)

        if returned:
            return self.meta.accept

    def download_data(self, client, new_sr=None, t1=None, t2=None,
                      returned=False, verbose=False):
        """
        Downloads seismograms based on a time interval.

        Parameters
        ----------
        client : :class:`~obspy.client.fdsn.Client`
            Client object
        new_sr : float
            New sampling rate (Hz)
        t1 : float
            Start time to consider (sec). Can be float or
            :class:`~obspy.core.UTCDateTime` object.
        t2 : float
            End time to consider (sec). Can be float or
            :class:`~obspy.core.UTCDateTime` object.
        returned : bool
            Whether or not to return the ``accept`` attribute
        verbose : bool
            Whether or not to print messages to screen during run-time


        Returns
        -------
        accept : bool
            Whether or not the object is accepted for further analysis

        Attributes
        ----------
        data : :class:`~obspy.core.Stream`
            Stream containing :class:`~obspy.core.Trace` objects

        """

        if self.meta is None:
            raise Exception("Requires event data as attribute - aborting")

        if not self.meta.accept:
            return

        try:
            # Define start time for request
            tstart = self.meta.time + t1
            tend = self.meta.time + t2

        except Exception:
            raise Exception(" Start and end request times have not been set")
            self.meta.accept = False
            if returned:
                return self.meta.accept

        # Get waveforms
        if verbose:
            print("* Requesting Waveforms: ")
            print("*    Startime: " + tstart.strftime("%Y-%m-%d %H:%M:%S"))
            print("*    Endtime:  " + tend.strftime("%Y-%m-%d %H:%M:%S"))

        # Download data
        err, stream = io.download_data(
            client=client, sta=self.sta, start=tstart, end=tend,
            zcomp=self.zcomp, verbose=verbose)

        if err or stream is None:
            print("*  Waveform Request Failed...Skipping")
            self.meta.accept = False
            if returned:
                return self.meta.accept

        # Store as attributes with traces in dictionary
        try:
            trE = stream.select(component='E')[0]
            trN = stream.select(component='N')[0]
            trZ = stream.select(component=self.zcomp)[0]
            self.data = Stream(traces=[trZ, trN, trE])

            # Filter Traces and resample
            if new_sr is not None:
                self.data.filter('lowpass', freq=0.5*new_sr,
                                 corners=2, zerophase=True)
                self.data.resample(new_sr)

        # If there is no ZNE, perhaps there is Z12?
        except Exception as e:
            tr1 = stream.select(component='1')[0]
            tr2 = stream.select(component='2')[0]
            trZ = stream.select(component=self.zcomp)[0]
            self.data = Stream(traces=[trZ, tr1, tr2])

            # Filter Traces and resample
            if new_sr is not None:
                self.data.filter('lowpass', freq=0.5*new_sr,
                                 corners=2, zerophase=True)
                self.data.resample(new_sr)

            if not np.sum([np.std(tr.data) for tr in self.data]) > 0.:
                print('Error: Data are all zeros')
                self.meta.accept = False

        if returned:
            return self.meta.accept

    def save(self, file):

        import pickle
        output = open(file, 'wb')
        pickle.dump(self, output)
        output.close()

        return


class BNG(Orient):
    """
    A BNG object inherits from :class:`~orientpy.classes.Orient`.
    This object contains a method to estimate station orientation
    based on P-wave particle motion.

    Notes
    -----

    This class is designed after the method developed by Braunmiller,
    Nabelek and Ghods (2020) [1]_. It is slightly more flexible, however, 
    in that it can handle either regional or teleseismic P-wave data
    and provides additional quality-control mesures (e.g., SNR, 1-R/Z).

    References
    ----------

    .. [1] Braunmiller, J., Nabelek, J., and Ghods, A., 2020, Sensor orientation 
       of Iranian broadband seismic stations from P-wave particle motion, 
       *Seismological Research Letters*, doi:10.1785/0220200019.

    Parameters
    ----------
    sta : :class:`~stdb.StDbElement`
        Object containing station information - from :mod:`~stdb` database.
    meta : :class:`~orientpy.classes.Meta`
        Object of metadata information for single event (initially set to None)
    data : :class:`~obspy.core.Stream`
        Stream object containing the three-component seismograms

    """

    def __init__(self, sta, zcomp='Z'):

        Orient.__init__(self, sta, zcomp=zcomp)

    def calc(self, dphi, dts, tt, bp=None, showplot=False):
        """
        Method to estimate azimuth of component `?H1` (or `?HN`). This
        method minimizes the energy (RMS) of the transverse component of
        P-wave data within some bandwidth.

        Parameters
        ----------
        dphi : float
            Azimuth increment for search (deg)
        dts : float
            Length of time window on either side of
            predicted P-wave arrival time (sec)
        tt : list
            List of two floats containing the time picks relative
            to P-wave time, within which to perform the estimation of
            station orientation (sec)
        bp : list
            List of two floats containing the low- and high-frequency
            corners of a bandpass filter (Hz)
        showplot : bool
            Whether or not to plot waveforms.

        Attributes
        ----------
        meta.phi : float
            Azimuth of H1 (or HN) component (deg)
        meta.cc : float
            Cross-correlation coefficient between
            vertical and radial component
        meta.snr : float
            Signal-to-noise ratio of P-wave measured on the
            vertical seismogram
        meta.TR : float
            Measure of the transverse to radial ratio. In reality
            this is 1 - T/R
        meta.RZ : float
            Measure of the radial to vertical ratio. In reality
            this is 1 - R/Z

        """

        # Work on a copy of the waveform data
        stream = self.data.copy()

        # Identify components in Stream
        comps_id = [tr.stats.channel[2] for tr in stream]
        test_sets = {
                        'ZNE': {'Z', 'N', 'E'},
                        'Z12': {'Z', '1', '2'},
                        'Z23': {'Z', '2', '3'},
                        '312': {'3', '1', '2'},
                        '123': {'1', '2', '3'}
                    }
        # Probably should raise an exception if set is '123',
        # as no correction is estimated for the vertical component

        for test_key in test_sets:
            test_set = test_sets[test_key]

            if test_set.issubset(set(comps_id)):

                # Use sets to avoid sorting complications
                comps_codes = list(test_key)
                break

        # Temporarily modify channel codes,
        # assuming that N/E are not oriented properly
        channel_code_prefix = stream[0].stats.channel[:2]

        stream.select(component=comps_codes[1])[0].stats.channel = channel_code_prefix + '1'
        stream.select(component=comps_codes[2])[0].stats.channel = channel_code_prefix + '2'
        stream.select(component=comps_codes[0])[0].stats.channel = channel_code_prefix + self.zcomp.upper()

        # Filter if specified
        if bp is not None:
            stream.filter('bandpass', freqmin=bp[0],
                          freqmax=bp[1], zerophase=True)

        # Get data and noise based on symmetric waveform wrt arrival
        start = stream[0].stats.starttime
        stnoise = stream.copy().trim(start, start+dts+tt[0])
        stdata = stream.copy().trim(start+dts+tt[0], start+dts+tt[1])

        # Define signal and noise
        tr1 = stdata.select(component='1')[0].copy()
        tr2 = stdata.select(component='2')[0].copy()
        trZ = stdata.select(component=self.zcomp.upper())[0].copy()
        ntrZ = stnoise.select(component=self.zcomp.upper())[0].copy()

        # Calculate and store SNR as attribute
        self.meta.snr = 10.*np.log10(
            utils.rms(trZ)*utils.rms(trZ)/utils.rms(ntrZ)/utils.rms(ntrZ))

        # Search through azimuths from 0 to 180 deg and find best-fit azimuth
        ang = np.arange(0., 180., dphi)
        cc1 = np.zeros(len(ang))
        cc2 = np.zeros(len(ang))
        cc3 = np.zeros(len(ang))
        cc4 = np.zeros(len(ang))
        for k, a in enumerate(ang):
            R, T = rotate_ne_rt(tr1.data, tr2.data, a)
            covmat = np.corrcoef(R, trZ.data)
            cc1[k] = covmat[0, 1]
            cc2[k] = 1. - utils.rms(T)/utils.rms(R)
            cc3[k] = utils.rms(T)
            cc4[k] = 1. - utils.rms(R)/utils.rms(trZ.data)

        # Get argument of minimum of cc3 and store useful measures
        ia = cc3.argmin()
        self.meta.cc = cc1[ia]
        self.meta.TR = cc2[ia]
        self.meta.RZ = cc4[ia]

        # correct for angles above 360
        phi = (self.meta.baz - float(ia)*dphi)

        # Use azimuth where CC is negative
        if self.meta.cc < 0.:
            phi += 180.

        if phi < 0.:
            phi += 360.
        if phi >= 360.:
            phi -= 360.

        # Store the best-fit azimuth
        self.meta.phi = phi

        # If a plot is requested, rotate Z12 to ZNE and then ZRT
        if showplot:

            # Now rotate components to proper N, E and then R, T
            sttmp = stream.copy()

            # Apply filter if defined previously
            if bp:
                sttmp.filter('bandpass', freqmin=bp[0],
                             freqmax=bp[1], zerophase=True)

            # Copy traces
            trN = sttmp.select(component='1')[0].copy()
            trE = sttmp.select(component='2')[0].copy()

            # Rotating from 1,2 to N,E is the negative of
            # rotation from RT to NE, with baz corresponding
            # to azim of component 1, or phi previously determined
            azim = self.meta.phi
            N, E = rotate_rt_ne(trN.data, trE.data, azim)
            trN.data = -1.*N
            trE.data = -1.*E

            # Update stats of streams
            trN.stats.channel = trN.stats.channel[:-1] + 'N'
            trE.stats.channel = trE.stats.channel[:-1] + 'E'

            # Store corrected traces in new stream and rotate to
            # R, T using back-azimuth
            stcorr = Stream(traces=[trN, trE])
            stcorr.rotate('NE->RT', back_azimuth=self.meta.baz)

            # Merge original and corrected streams
            st = stream + stcorr

            # Plot
            plot = plotting.plot_bng_waveforms(self, st, dts, tt)
            plot.show()

        return


class DL(Orient):
    """
    A DL object inherits from :class:`~orientpy.classes.Orient`. This object
    contains a method to estimate station orientation based on Rayleigh-wave
    particle motion.

    Notes
    -----

    This algorithm is heavily based on the code `DLOPy` by
    Doran and Laske (2017) [2]_. The original code can be found here:
    https://igppweb.ucsd.edu/~adoran/DLOPy.html

    There are also a couple of versions available on GitHub:

    - https://github.com/kschramm-usgs/DLOPy
    - https://github.com/jbrussell/DLOPy_v1.0

    References
    ----------

    .. [2] Doran, A. K., and G. Laske (2017). Ocean-bottom seismometer
       instrument orientations via automated Rayleigh-wave arrival-angle
       measurements, *Bulletin of the Seismological Society of America*,
       107(2), doi:10.1785/0120160165.

    Parameters
    ----------
    sta : :class:`~stdb.StDbElement`
        Object containing station information - from :mod:`~stdb` database.
    meta : :class:`~orientpy.classes.Meta`
        Object of metadata information for single event (initially set to None)
    data : :class:`~obspy.core.Stream`
        Stream object containing the three-component seismograms

    """

    def __init__(self, sta, zcomp='Z'):

        Orient.__init__(self, sta, zcomp=zcomp)

    def calc(self, showplot=False):
        """
        Method to estimate azimuth of component `?H1` (or `?HN`). This
        method maximizes the normalized covariance between the
        Hilbert-transformed vertical component and the radial component of
        Rayleigh-wave data measured at multiple periods. This is done for the
        shortest (R1) and longest (R2) orbits of fundamental-mode Rayleigh
        waves. Window selection is done based on average group velocity
        extracted from a global model of Rayleigh-wave dispersion.

        Parameters
        ----------
        showplot : bool
            Whether or not to plot waveforms.

        Attributes
        ----------
        meta.R1phi : list
            List of azimuth of H1 (or HN) component (deg) based on R1, measured
            at different periods
        meta.R1cc : float
            Cross-correlation coefficient between hilbert-transformed vertical
            and radial component for R1, measured at different periods
        meta.R2phi : list
            List of azimuth of H1 (or HN) component (deg) based on R2, measured
            at different periods
        meta.R2cc : float
            Cross-correlation coefficient between hilbert-transformed vertical
            and radial component for R2, measured at different periods

        """

        # Work on a copy of the waveform data
        stream = self.data.copy()

        # Initialize surface wave arrays
        nper = 7
        R1phi = np.zeros(nper)
        R1cc = np.zeros(nper)
        R2phi = np.zeros(nper)
        R2cc = np.zeros(nper)

        # Load group velocity maps
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

        # Get parameters for R2
        Rearth = 6371.25
        circE = 2.*np.pi*Rearth
        dist2 = circE - self.meta.epi_dist
        baz2 = self.meta.baz + 180.
        if baz2 >= 360.:
            baz2 -= 360.

        # Check data length, data quality
        if utils.checklen(self.data, 4.*60.*60.):
            raise Exception("      Error: Length Incorrect")

        # Get path-averaged group velocities
        Ray1, Ray2 = utils.pathvels(
            self.sta.latitude, self.sta.longitude,
            self.meta.lat, self.meta.lon,
            map10, map15, map20, map25, map30, map35, map40)

        # Calculate arrival angle for each frequency and orbit
        Rf = [40., 35., 30., 25., 20., 15., 10.]
        LPF = [0.035, 0.03, 0.025, 0.02, 0.015, 0.01, 0.005]
        HPF = [0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015]
        winlen1 = [20., 17., 14., 12., 10., 10., 7.]
        winlen2 = [24., 20., 16., 13., 10., 10., 7.]
        flist = list(zip(Rf, HPF, LPF, winlen1, winlen2))
        for k, item in enumerate(flist):

            # R1 path
            R1phi[k], R1cc[k] = utils.DLcalc(
                stream, item[0], item[1],
                item[2], self.meta.epi_dist, self.meta.baz, Ray1,
                winlen=item[3], ptype=0, zcomp=self.zcomp)

            # R2 path
            R2phi[k], R2cc[k] = utils.DLcalc(
                stream, item[0], item[1],
                item[2], dist2, baz2, Ray2,
                winlen=item[4], ptype=0, zcomp=self.zcomp)

        # Store azimuths and CC values as attributes
        self.meta.R1phi = R1phi
        self.meta.R2phi = R2phi
        self.meta.R1cc = R1cc
        self.meta.R2cc = R2cc

        return
