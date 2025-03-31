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
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from orientpy import utils
import numpy as np
from scipy.stats import gaussian_kde

def density_estimate(values, x, alpha):
    """ 
    This function estimates the KDE, mode of distribution, and
    95% confidence intervals for the ``values`` variables
    evaluated at locations ``x``

    Parameters
    ----------
    values : :class:`~numpy.ndarray`
        Array of values for which to estimate KDE
    x : :class:`~numpy.ndarray`
        Array of locations at which to evaluate the KDE
    alpha : float
        Significance level (alpha=0.05 means 95% confidence)

    Returns
    -------
    kde : :class:`~scipy.stats.gaussian_kde`
        KDE object for data ``values`` evaluated at ``x``
    x_max : float
        Value at which the KDE is maximum
    CI_min : float
        Lower bound of confidence interval
    CI_max : float
        Upper bound of confidence interval
    
    """
    alpha /= 2.
    kernel = gaussian_kde(values)
    kde = kernel.evaluate(x)
    dx = x[1] - x[0]
    cdf = np.cumsum(kernel.evaluate(x)*dx)
    x_max = np.mean(x[kernel(x) == np.max(kernel(x))])
    CI_min = x[cdf < alpha][-1]
    try:
        CI_max = x[cdf > (1.-alpha)][0]
    except:
        CI_max = np.max(x)

    return kde, x_max, CI_min, CI_max


def plot_bng_waveforms(bng, stream, dts, tt):
    """ 
    This function plots the original and rotated waveforms
    following the BNG processing for quality control.

    Parameters
    ----------
    bng : :class:`~orientpy.BNG`
        BNG object for a single event
    stream : :class:`~obspy.core.Stream`
        Stream of un-rotated and rotated waveforms
    dts : float
        Length of time window over which to plot waveforms
    tt : list
        Zoomed in time window over which processing is done

    Returns
    -------
    plt : :class:`~matplotlib.pyplot`
        Handle to final plot

    """

    fig, ax = plt.subplots(5, 1, sharey=True)
    cmpts = [bng.zcomp.upper(), 'R', 'T', '1', '2']
    taxis = np.linspace(-dts, dts, stream[0].stats.npts)

    for item in list(zip(ax, cmpts)):
        item[0].plot(taxis, stream.select(
            component=item[1])[0].data, 'k', lw=1.)
        item[0].axvspan(tt[0], tt[1], color='coral', alpha=0.5)
        item[0].axvline(0., ls='--', c='k')
        item[0].set_ylabel('H'+item[1])

    ax[-1].set_xlabel('Time following predicted P-wave arrival (sec)')

    text = "Station "+bng.sta.station+": " + \
        "{0:.1f} deg\n".format(bng.meta.phi) + \
        "SNR: {0:.1f}, CC: {1:.1f}, ".format(
            bng.meta.snr, np.abs(bng.meta.cc)) + \
        "1-T/R: {0:.1f}, 1-R/Z: {1:.1f}".format(
            bng.meta.TR, bng.meta.RZ)
    plt.suptitle(text, fontsize=10)
    plt.subplots_adjust(hspace=0.05)

    return plt


def plot_bng_conditions(stkey, snr, cc, TR, RZ, ind):
    """ 
    This function plots all parameters and the threshold condition for
    contributing to the final estimate.

    Parameters
    ----------
    stkey : str
        Station key
    snr : :class:`~numpy.ndarray`
        Array of signal-to-noise ratio values for each earthquake
    cc : :class:`~numpy.ndarray`
        Array of cross-correlation values between rotated radial 
        and vertical components
    TR : :class:`~numpy.ndarray`
        Array of transverse-to-radial component ratios for rotated
        components
    RZ : :class:`~numpy.ndarray`
        Array of radial-to-vertical component ratios for rotated
        components
    ind : :class:`~numpy.ndarray`
        Array of boolean (index) values where conditions on previous 
        parameters are examined

    Returns
    -------
    plt : :class:`~matplotlib.pyplot`
        Handle to final plot

    """

    f = plt.figure(figsize=(7.5, 5))
    gs = gridspec.GridSpec(2, 3)
    gs.update(wspace=0.45, hspace=0.3)

    # Ax1: phi vs snr
    ax1 = f.add_subplot(gs[0])
    ax1.scatter(snr, cc, marker='+')
    ax1.scatter(snr[ind], cc[ind], marker='+')
    ax1.set_xlabel('SNR', fontsize=10)
    ax1.set_ylabel('CC', fontsize=10)
    ax1.tick_params(axis='both', labelsize=10)

    # Ax2: phi vs CC
    ax2 = f.add_subplot(gs[1])
    ax2.scatter(snr, TR, marker='+')
    ax2.scatter(snr[ind], TR[ind], marker='+')
    ax2.set_xlabel('SNR', fontsize=10)
    ax2.set_ylabel('1 - T/R', fontsize=10)

    # Ax3: phi vs TR
    ax3 = f.add_subplot(gs[2])
    ax3.scatter(snr, RZ, marker='+')
    ax3.scatter(snr[ind], RZ[ind], marker='+')
    ax3.set_xlabel('SNR', fontsize=10)
    ax3.set_ylabel('1 - R/Z', fontsize=10)

    # Ax4: phi vs RZ
    ax4 = f.add_subplot(gs[3])
    ax4.scatter(cc, TR, marker='+')
    ax4.scatter(cc[ind], TR[ind], marker='+')
    ax4.set_xlabel('CC', fontsize=10)
    ax4.set_ylabel('1 - T/R', fontsize=10)

    # Ax5: phi vs baz
    ax5 = f.add_subplot(gs[4])
    ax5.scatter(cc, RZ, marker='+')
    ax5.scatter(cc[ind], RZ[ind], marker='+')
    ax5.set_xlabel('CC', fontsize=10)
    ax5.set_ylabel('1 - R/Z', fontsize=10)

    # Ax6: phi vs mag
    ax6 = f.add_subplot(gs[5])
    ax6.scatter(TR, RZ, marker='+')
    ax6.scatter(TR[ind], RZ[ind], marker='+')
    ax6.set_xlabel('1 - T/R', fontsize=10)
    ax6.set_ylabel('1 - R/Z', fontsize=10)

    text = "Station "+stkey
    plt.suptitle(text, fontsize=12)
    # gs.tight_layout(f, rect=[0, 0, 1, 0.95])

    return plt


def plot_bng_results(stkey, phi, snr, cc, TR, RZ, baz, mag,
                 ind, val, err, alpha=0.05):
    """ 
    This function plots the results of all BNG estimates with final
    estimates from those that pass the conditions.

    Parameters
    ----------
    stkey : str
        Station key
    phi : :class:`~numpy.ndarray`
        Array of estimated azimuth for each earthquake
    snr : :class:`~numpy.ndarray`
        Array of signal-to-noise ratio values for each earthquake
    cc : :class:`~numpy.ndarray`
        Array of cross-correlation values between rotated radial 
        and vertical components
    TR : :class:`~numpy.ndarray`
        Array of transverse-to-radial component ratios for rotated
        components
    RZ : :class:`~numpy.ndarray`
        Array of radial-to-vertical component ratios for rotated
        components
    baz : :class:`~numpy.ndarray`
        Array of back-azimuth values corresponding to each earthquake
    mag : :class:`~numpy.ndarray`
        Array of magnitude values corresponding to each earthquake
    ind : :class:`~numpy.ndarray`
        Array of boolean values where conditions on previous parameters
        are examined
    val : float
        Final estimated azimuth for station orientation
    err : float
        Final error estimate on azimuth
    alpha : float
        Significance level (alpha=0.05 means 95% confidence)

    Returns
    -------
    plt : :class:`~matplotlib.pyplot`
        Handle to final plot

    """

    # Re-center phi values and extract ones that meet condition
    allphi = utils.centerat(phi, m=val)
    goodphi = utils.centerat(phi[ind], m=val)

    # Figure
    fig = plt.figure(figsize=(8, 5))

    gs = gridspec.GridSpec(1, 2, width_ratios=[8, 1])
    gs.update(wspace=0.1)
    gs0 = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=gs[0])
    gs1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1])

    data = [snr, cc, TR, RZ, baz, mag]
    xlab = ['SNR', 'CC', '1 - T/R', '1 - R/Z',
            'Back-azimuth ($^\circ$)', 'Magnitude']

    for item in list(zip(gs0, data, xlab)):

        x = np.linspace(item[1].min(), item[1].max(), 20)
        ax = fig.add_subplot(item[0])
        ax.fill_between(x, val+err, val-err, fc='grey', alpha=0.5)
        ax.axhline(val, lw=1, ls='-', c='k')
        ax.scatter(item[1], allphi, marker='+')
        ax.scatter(item[1][ind], goodphi, marker='+')
        ax.set_xlabel(item[2], fontsize=10)
        ax.tick_params(axis='both', labelsize=10)
        ax.set_ylim([val-180., val+180.])

    fig.axes[0].set_ylabel('BH1 Orientation \n Angle ($^\circ$)', fontsize=10)
    fig.axes[3].set_ylabel('BH1 Orientation \n Angle ($^\circ$)', fontsize=10)

    text = "Station "+stkey + \
        ": $\phi$ = {0:.1f} $\pm$ {1:.1f}".format(val, err)
    plt.suptitle(text, fontsize=12)

    # KDE plot
    y = np.arange(val-180., val+180., 0.01)

    phip = utils.outlier(allphi, 5.)
    phipp = utils.outlier(goodphi, 5.)

    if len(phip) > 0 and len(phipp) > 0:

        k1, k1_map, CI1_min, CI1_max = density_estimate(phip, y, alpha)
        k2, k2_map, CI2_min, CI2_max = density_estimate(phipp, y, alpha)

        ax = fig.add_subplot(gs1[0], sharey=fig.axes[2])
        ax.plot(k1, y)
        ax.plot(k2, y)
        ax.yaxis.tick_right()

        ax = fig.add_subplot(gs1[1], sharey=fig.axes[-1])
        ax.plot(k1, y)
        ax.plot(k2, y)
        ax.yaxis.tick_right()

        ax.set_xlabel('KDE')

    # plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.subplots_adjust(hspace=0.27, wspace=0.3)

    return plt


def plot_dl_results(stkey, R1phi, R1cc, R2phi, R2cc, ind, val,
                 err, phi, cc, cc0, alpha=0.05):
    """ 
    This function plots the results of all DL estimates with final
    estimates from those that pass the conditions.

    Parameters
    ----------
    stkey : str
        Station key
    R1phi : :class:`~numpy.ndarray`
        Array of azimuth values from each estimate for direct 
        Rayleigh-wave pass
    R1cc : :class:`~numpy.ndarray`
        Array of cross-correlation values between rotated radial 
        and vertical components for direct Rayleigh-wave pass
    R2phi : :class:`~numpy.ndarray`
        Array of azimuth values from each estimate for complementary 
        Rayleigh-wave pass
    R2cc : :class:`~numpy.ndarray`
        Array of cross-correlation values between rotated radial 
        and vertical components for complementary Rayleigh-wave pass
    ind : :class:`~numpy.ndarray`
        Array of boolean values where conditions on previous parameters
        are examined
    val : float
        Final estimated azimuth for station orientation
    err : float
        Final error estimate on azimuth
    phi : :class:`~numpy.ndarray`
        All azimuth estimates from both R1 (direct pass) and R2 
        (complementary pass)
    cc : :class:`~numpy.ndarray`
        All cross-correlation estimates from both R1 (direct pass) and 
        R2 (complementary pass)
    cc0 : float
        Threshold cross-correlation value
    alpha : float
        Significance level (alpha=0.05 means 95% confidence)

    Returns
    -------
    plt : :class:`~matplotlib.pyplot`
        Handle to final plot

    """

    # Re-center phi values and extract ones that meet condition
    allphi = utils.centerat(phi, m=val)
    goodphi = utils.centerat(phi[ind], m=val)

    f = plt.figure(figsize=(8, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])

    # Scatter plot
    x = np.arange(0., 1.01, 0.01)

    ax1 = f.add_subplot(gs[0])
    ax1.fill_between(x, val+err, val-err,
                     fc='grey', alpha=0.5)
    ax1.axvline(cc0, c='k', ls='--', lw=1.)
    ax1.axhline(val, lw=1, ls='-', c='k')
    ax1.scatter(R1cc, utils.centerat(R1phi, m=val),
                marker='x', label='R1')
    ax1.scatter(R2cc, utils.centerat(R2phi, m=val),
                marker='+', label='R2')
    ax1.set_ylabel('BH1 Orientation \n Angle ($^\circ$)', fontsize=10)
    ax1.set_xlabel('Cross-correlation value', fontsize=10)
    ax1.set_ylim([val-180., val+180.])
    ax1.set_xlim([0, 1])
    ax1.tick_params(axis='both', labelsize=10)
    ax1.legend()

    # KDE plot
    y = np.arange(val-180., val+180., 0.01)

    # Apply outlier analysis
    phip = utils.outlier(allphi, 5.)
    phipp = utils.outlier(goodphi, 5.)

    if len(phip) > 0 and len(phipp) > 0:

        k1, k1_map, CI1_min, CI1_max = density_estimate(phip, y, alpha)
        k2, k2_map, CI2_min, CI2_max = density_estimate(phipp, y, alpha)

        ax2 = f.add_subplot(gs[1], sharey=ax1)
        ax2.plot(k1, y, label='All')
        ax2.plot(k2, y, label='Robust')

        ax2.set_xlabel('KDE')
        ax2.yaxis.tick_right()
        ax2.legend()

    # Add text as title
    text = "Station "+stkey + \
        ": $\phi$ = {0:.1f} $\pm$ {1:.1f}  (n={2:.0f}, cc={3:.2f})".format(val, err, sum(ind),cc0)
    plt.suptitle(text, fontsize=12)
    gs.tight_layout(f, rect=[0, 0, 1, 0.95])

    return plt
    
