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
import numpy as np
import matplotlib.pyplot as plt
from geographiclib.geodesic import *
from obspy.clients.fdsn import Client
from obspy.signal.rotate import rotate_rt_ne, rotate_ne_rt
from obspy import UTCDateTime, read, Stream, Trace
from scipy.stats import circmean as cmean
from scipy.stats import circstd as cstd
from scipy.stats import hmean as hm
import scipy.signal as sig
import os
import math


def catclean(cat):
    """ 
    This function looks for repeat events in a catalogue of 
    earthquakes.

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    cat : :class:`~obspy.core.event.Catalog`
        Catalogue XML object, container for Event objects.

    Returns
    -------
    reps : :class:`~numpy.ndarray`
        Array of ID of repeat events in catalogue.
    
    """

    def close(x1, x2, val):
        if np.abs(x1 - x2) < val:
            return True
        else:
            return False

    rep = np.array([], dtype=int)
    for k, event in enumerate(cat):
        for kk, ev in enumerate(cat):
            if (ev.origins[0].time - event.origins[0].time < 60.*60.) and \
                    close(ev.origins[0].latitude,
                          event.origins[0].latitude, 0.8) and \
                    close(ev.origins[0].longitude,
                          event.origins[0].longitude, 0.8) and \
                    (ev.magnitudes[0].mag < event.magnitudes[0].mag - 0.3):
                rep = np.append(rep, kk)
    return rep


def checklen(st, hrs):
    """
    Function to check if there is enough downloaded data to run program

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    st : :class:`~obspy.core.Stream`
        Stream containing waveforms
    hrs : float
        Amount of time (hours) that should be available for analysis

    Returns
    -------

    boolean
        Whether or not the check is successful

    """
    L = len(st)
    for i in np.arange((L)):
        if (UTCDateTime(st[i].stats.endtime) -
                UTCDateTime(st[i].stats.starttime)) + 100. < hrs:
            return True
    if np.var(st[0].data) < 1 or np.var(st[1].data) < 1 or \
            np.var(st[2].data) < 1:
        return True
    return False


def resiz(x1, x2, x3):
    """
    Function to resize arrays to all identical shapes

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    x* : :class:`~numpy.ndarray`
        Array for which to check length (here it's trace.data)

    Returns
    -------

    x* : :class:`~numpy.ndarray`
        Array trimmed to shortest segment

    """
    L = np.min(np.array([len(x1), len(x2), len(x3)]))
    return x1[0:L], x2[0:L], x3[0:L]


def getf(freq, A):
    """
    Function to extract frequency in array.

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    freq : float
        Frequency of interest (Hz)
    A : :class:`~numpy.ndarray`
        Array of frequencies obtained from global dispersion maps

    Returns
    -------
    v : float
        Nearest frequency in A

    """

    for i in np.arange((len(A))):
        if A[i, 0] == freq:
            v = A[i, 1]
    return v


def nv(x, v):
    """
    Function to find nearest value in an np.array.

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    x : :class:`~numpy.ndarray`
        Input array
    v : float
        Value to compare within array

    Returns
    -------
    x[idx] : float
        Nearest value to v in x

    """

    return x[(np.abs(x - v)).argmin()]


def rms(x):
    """
    Function to calculate root-mean-square of array

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    x : :class:`~numpy.ndarray`
        Input array

    Returns
    -------
    rms : float
        Root-Mean-Square value of `x`

    """

    return np.sqrt(np.mean(np.abs(x)**2))


def mad(x):
    """
    Function to calculate Median Absolute Deviation

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    x : :class:`~numpy.ndarray`
        Input array

    Returns
    -------
    mad : float
        Median Absolute deviation of `x`

    """

    return np.median(np.abs(x - np.median(x)))


def outlier(x, lim=5.0):
    """
    Function to remove outliers based on MAD threshold

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    x : :class:`~numpy.ndarray`
        Input array

    Returns
    -------
    x : :class:`~numpy.ndarray`
        Shortened array where MAD is lower than `lim`

    """

    return x[np.abs(x - np.median(x))/mad(x) < lim]


def boot(x, bootnum):
    """
    Function to calculate directional mean value from bootstrap

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    x : :class:`~numpy.ndarray`
        Input array
    bootnum : int
        Number of bootstrap estimates to produce

    Returns
    -------
    m : :class:`~numpy.ndarray`
        Array of directional mean values

    """

    m = np.zeros((bootnum))
    L = len(x)
    for i in np.arange((bootnum)):
        a = np.random.choice(x, size=L, replace=True)
        m[i] = cmean(a, high=360)
    return m


def centerat(phi, m=0.):
    """
    Function to re-center azimuth array to the mean value

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    phi : :class:`~numpy.ndarray`
        Input array
    m : float
        Value around which to center the azimuth

    Returns
    -------
    phinew : :class:`~numpy.ndarray`
        Re-centered azimuth array

    """

    phinew = np.copy(phi)
    if len(phi.shape) == 1:
        for i in np.arange((len(phi))):
            if phi[i] >= m+180.:
                phinew[i] -= 360.
            if phi[i] <= m-180.:
                phinew[i] += 360.
    else:
        for k in np.arange((phi.shape[1])):
            for i in np.arange((phi.shape[0])):
                if phi[i, k] >= m+180.:
                    phinew[i, k] -= 360.
                if phi[i, k] <= m-180.:
                    phinew[i, k] += 360.
    return phinew


def pathvels(lat1, lon1, lat2, lon2,
             map10, map15, map20, map25, map30, map35, map40):
    """
    Overall function to get path-averaged group velocity.

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    lat1 : float
        Latitude of origin point (deg)
    lon1 : float
        Longitude of origin point (deg)
    lat2 : float
        Latitude of end point (deg)
    lon2 : float
        Longitude of end point (deg)
    map* : :class:`~numpy.ndarray`
        maps of Rayleigh-wave dispersion at various frequencies

    Returns
    -------
    R1 : :class:`~numpy.ndarray`
        R1 velocity path 

    R2 : :class:`~numpy.ndarray`
        R2 velocity path 

    """

    Rearth = 6371.25
    circE = 2.*np.pi*Rearth

    # Get distance and azimuth
    p = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)

    minor = float(p['s12']) / 1000.
    major = circE-minor

    l1 = Geodesic.WGS84.Line(lat1, lon1, p['azi1'])
    l2 = Geodesic.WGS84.Line(lat1, lon1, p['azi1'] - 180.)

    deg = 111.17
    D1 = np.zeros([1, 2])
    for i in np.arange((361)):
        b = l1.Position(deg*1000*i)
        p2 = np.array([b['lat2'], b['lon2']])

        if i == 0:
            D1[0, 0] = p2[0]
            D1[0, 1] = p2[1]
        else:
            D1 = np.vstack((D1, p2))

        bb = Geodesic.WGS84.Inverse(lat2, lon2, b['lat2'], b['lon2'])
        if bb['s12'] <= deg*1000.:
            break

    D2 = np.zeros([1, 2])
    for i in np.arange((361)):
        b = l2.Position(deg*1000*i)
        p2 = np.array([b['lat2'], b['lon2']])

        if i == 0:
            D2[0, 0] = p2[0]
            D2[0, 1] = p2[1]
        else:
            D2 = np.vstack((D2, p2))

        bb = Geodesic.WGS84.Inverse(lat2, lon2, b['lat2'], b['lon2'])
        if bb['s12'] <= deg*1000.0:
            break

    # We now have lat and lon points along the major and minor great circles.
    # We calculate the group velocity of at each point, and then
    # find the average velocities.
    for k in np.arange((len(D1))):
        if D1[k, 1] < 0:
            D1[k, 1] += 360
    for k in np.arange((len(D2))):
        if D2[k, 1] < 0:
            D2[k, 1] += 360

    def Ray(D):
        U1 = np.zeros([len(D), 7])
        for k in np.arange((len(D))):
            # Do latitude first

            # Get correct precision
            # designed to match results of Ma et al codes
            if np.abs(D[k, 1]) < 10:
                D[k, 1] = round(D[k, 1], 5)
            if np.abs(D[k, 1]) >= 10 and np.abs(D[k, 1]) < 100:
                D[k, 1] = round(D[k, 1], 4)
            if np.abs(D[k, 1]) > 100:
                D[k, 1] = round(D[k, 1], 3)
            if np.abs(D[k, 0]) < 10:
                D[k, 0] = round(D[k, 0], 5)
            if np.abs(D[k, 0]) >= 10:
                D[k, 0] = round(D[k, 0], 4)

            # find right index
            q = np.where(map10[:, 1] == nv(map10[:, 1], (D[k, 0])))[0]
            qq = np.where(map10[q, 0] == nv(map10[q, 0],  (D[k, 1])))[0]
            idx = q[qq]

            # update path
            U1[k, 0] = map10[idx, 2]
            U1[k, 1] = map15[idx, 2]
            U1[k, 2] = map20[idx, 2]
            U1[k, 3] = map25[idx, 2]
            U1[k, 4] = map30[idx, 2]
            U1[k, 5] = map35[idx, 2]
            U1[k, 6] = map40[idx, 2]
        mhz = np.array([10, 15, 20, 25, 30, 35, 40])
        return np.array((mhz, hm(U1, axis=0))).T

    return Ray(D1), Ray(D2)


def DLcalc(stream, Rf, LPF, HPF, epi, baz, A, winlen=10., ptype=0):
    """
    DORAN-LASKE calculation for one freq, one orbit of surface wave

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    stream : float
        Latitude of origin point (deg)
    lon1 : float
        Longitude of origin point (deg)
    lat2 : float
        Latitude of end point (deg)
    lon2 : float
        Longitude of end point (deg)
    map* : :class:`~numpy.ndarray`
        maps of Rayleigh-wave dispersion at various frequencies

    Returns
    -------
    R1 : :class:`~numpy.ndarray`
        R1 velocity path 

    R2 : :class:`~numpy.ndarray`
        R2 velocity path 

    """

    # Pre-process
    stream.taper(type='hann', max_percentage=0.05)
    stream.filter("lowpass", freq=LPF, corners=4, zerophase=True)
    stream.filter("highpass", freq=HPF, corners=4, zerophase=True)
    stream.detrend()

    # Window info
    Rvel = getf(Rf, A)  # Group velocity at Rf
    R1window = (1.0/(Rf/1000.))*winlen
    arv = 1./Rvel * epi
    r1 = arv - R1window/2.
    r2 = arv + R1window/2.

    dt = stream[0].stats.starttime
    st = stream.slice(starttime=dt+r1, endtime=dt+r2)

    # Extract waveform data for each component
    try:
        tr1 = st.select(component='1')[0].data
        tr2 = st.select(component='2')[0].data
    except:
        tr1 = st.select(component='N')[0].data
        tr2 = st.select(component='E')[0].data
    trZ = st.select(component='Z')[0].data

    # Calculate Hilbert transform of vertical trace data
    trZ = np.imag(sig.hilbert(trZ))

    # Ensure all data vectors are same length
    tr1, tr2, trZ = resiz(tr1, tr2, trZ)

    # Rotate through and find max normalized covariance
    dphi = 0.1
    ang = np.arange(0., 360., dphi)
    cc1 = np.zeros(len(ang))
    cc2 = np.zeros(len(ang))
    for k, a in enumerate(ang):
        R, T = rotate_ne_rt(tr1, tr2, a)
        covmat = np.corrcoef(R, trZ)
        cc1[k] = covmat[0, 1]
        cstar = np.cov(trZ, R)/np.cov(trZ)
        cc2[k] = cstar[0, 1]

    # Get argument of maximum of cc2
    ia = cc2.argmax()  

    # Get azimuth and correct for angles above 360
    phi = (baz - float(ia)*dphi) + 180.
    if phi < 0.:
        phi += 360.
    if phi >= 360.:
        phi -= 360.

    # # plotting:
    # # ptype=0, no plot
    # # ptype=1, Rayleigh plot
    # # ptype=2, Love plot
    # if ptype == 1:
    #     import matplotlib.dates as dat
    #     X = P[0].times()
    #     T = np.zeros((len(X)))
    #     for q in np.arange((len(T))):
    #         T[q] = dt + r1 + X[q]
    #     ZZ = dat.epoch2num(T)
    #     Z = dat.num2date(ZZ)
    #     n, e = rot2d(rdat, rdat2, ANG/4.)
    #     plt.figure()
    #     plt.plot(Z, vdat, label='Vertical')
    #     plt.hold("on")
    #     plt.plot(Z, n, label='BH1')
    #     plt.legend(loc=4)
    #     plt.xlabel('Time')
    #     plt.ylabel('Counts')
    #     plt.title('D-L Results (%i mHz)' % (Rf))
    # elif ptype == 2:
    #     import matplotlib.dates as dat
    #     X = P[0].times()
    #     T = np.zeros((len(X)))
    #     for q in np.arange((len(T))):
    #         T[q] = dt+r1+X[q]
    #     ZZ = dat.epoch2num(T)
    #     Z = dat.num2date(ZZ)
    #     n, e = rot2d(rdat, rdat2, ANG/4.)
    #     plt.figure()
    #     plt.subplot(121)
    #     plt.plot(Z, vdat, label='Vertical')
    #     plt.hold("on")
    #     plt.plot(Z, n, label='BH1')
    #     plt.legend(loc=4)
    #     plt.xlabel('Time')
    #     plt.suptitle('D-L Results (%i mHz)' % (Rf))
    #     plt.subplot(122)
    #     plt.plot(Z, e, label='BH2')
    #     plt.xlabel('Time')
    #     plt.ylabel('Counts')
    #     plt.legend(loc=4)
    # elif ptype == 3:
    #     import matplotlib.dates as dat
    #     X = P[0].times()
    #     T = np.zeros((len(X)))
    #     for q in np.arange((len(T))):
    #         T[q] = dt+r1+X[q]
    #     n, e = rot2d(rdat, rdat2, ANG/4.)
    #     plt.figure()
    #     plt.plot(T, vdat, label='Vertical')

    # if ptype > 0:
    #     plt.show()

    return phi, cc1[ia]


# Final calculation
def estimate(phi, ind):
    """
    Function to estimate final azimuth from 

    ADRIAN. K. DORAN and GABI LASKE, DLOPy VERSION 1.0, 
    RELEASED APRIL 2017

    Parameters
    ----------
    phi : :class:`~numpy.ndarray`
        Input array of estimated azimuth
    ind : list
        List of index values that satisfy some QC condition 
        for phi

    Returns
    -------
    m : float
        Mean value of robust, bootstrapped estimates

    s : float
        2 Standard deviations of robust, bootstrapped estimates

    """

    # Get phi', where conditions are met
    phip = phi[ind]

    # Re-center at circular mean
    phip = centerat(phip, m=cmean(phi, high=360))

    # Remove outliers
    phipp = outlier(phip, 5.)

    # bootstrap results for statistic
    m = boot(phipp, 5000)

    return cmean(m, high=360), 2*1.96*cstd(m, high=360)


