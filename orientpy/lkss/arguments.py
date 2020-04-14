# Copyright 2020 Pascal Audet
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
"""

Module containing the main utility functions used in the `RfPy` scripts
that accompany this package.

"""

# -*- coding: utf-8 -*-
from obspy import UTCDateTime
from numpy import nan, isnan
from obspy.core import Stream, read
from os.path import exists as exist
from argparse import ArgumentParser


def get_lkss_arguments():
    """
    Get Options from :class:`~argparse.ArgumentParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="Usage: %(prog)s [options] <station database>",
        description="Script used to find orientation of station from "+
        "receiver function data ")

    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)
    # General Settings
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys for " +
        "which to perform the analysis. These must be " +
        "contained within the station database. Partial keys will " +
        "be used to match against those in the dictionary. For " +
        "instance, providing IU will match with all stations in " +
        "the IU network [Default processes all stations in the database]")
    parser.add_argument(
        "-v", "-V", "--verbose",
        action="store_true",
        dest="verb",
        default=False,
        help="Specify to increase verbosity.")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing figures. " +
        "[Default False]")

    PreGroup = parser.add_argument_group(
        title='Pre-processing Settings',
        description="Options for pre-processing of receiver function " +
        "data before plotting")
    PreGroup.add_argument(
        "--snr",
        action="store",
        type=float,
        dest="snr",
        default=-9999.,
        help="Specify the SNR threshold for extracting receiver functions. " +
        "[Default None]")
    PreGroup.add_argument(
        "--snrh",
        action="store",
        type=float,
        dest="snrh",
        default=-9999.,
        help="Specify the horizontal component SNR threshold for extracting receiver functions. " +
        "[Default None]")
    PreGroup.add_argument(
        "--cc",
        action="store",
        type=float,
        dest="cc",
        default=-1.,
        help="Specify the CC threshold for extracting receiver functions. " +
        "[Default None]")
    PreGroup.add_argument(
        "--no-outlier",
        action="store_true",
        dest="no_outl",
        default=False,
        help="Set this option to delete outliers based on the MAD "+
        "on the variance. [Default False]")
    PreGroup.add_argument(
        "--bp",
        action="store",
        type=str,
        dest="bp",
        default=None,
        help="Specify the corner frequencies for the bandpass filter. " +
        "[Default no filtering]")
    PreGroup.add_argument(
        "--pws",
        action="store_true",
        dest="pws",
        default=False,
        help="Set this option to use phase-weighted stacking during binning "+
        " [Default False]")
    PreGroup.add_argument(
        "--nbaz",
        action="store",
        dest="nbaz",
        type=int,
        default=72,
        help="Specify integer number of back-azimuth bins to consider " +
        "(typically 36 or 72). If not None, the plot will show receiver " +
        "functions sorted by back-azimuth values. [Default 72]")
    PreGroup.add_argument(
        "--trange",
        action="store",
        default=None,
        type=str,
        dest="trange",
        help="Specify the time range for decomposition (sec). Negative times "+
        "are allowed [Default -1., 1.]")
    PreGroup.add_argument(
        "--boot",
        action="store_true",
        dest="boot",
        default=False,
        help="Set this option to calculate bootstrap statistics "+
        " [Default False]")
    PreGroup.add_argument(
        "--plot-f",
        action="store_true",
        dest="plot_f",
        default=False,
        help="Set this option to plot the function f(phi) "+
        "[Default False]")
    PreGroup.add_argument(
        "--plot-comps",
        action="store_true",
        dest="plot_comps",
        default=False,
        help="Set this option to plot the misoriented and rotated harmonic "+
        "components [Default False]")

    args = parser.parse_args()

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    if args.bp is not None:
        args.bp = [float(val) for val in args.bp.split(',')]
        args.bp = sorted(args.bp)
        if (len(args.bp)) != 2:
            parser.error(
                "Error: --bp should contain 2 " +
                "comma-separated floats")

    if args.trange is None:
        args.tmin = -1.
        args.tmax = 1.
    if args.trange is not None:
        args.trange = [float(val) for val in args.trange.split(',')]
        args.trange = sorted(args.trange)
        if (len(args.trange)) != 2:
            parser.error(
                "Error: --trange should contain 2 " +
                "comma-separated floats")

    return args

