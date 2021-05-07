#!/usr/bin/env python

# Final Orientation Calculation File

import stdb
import pickle
import os.path
import numpy as np
from orientpy import utils, plotting
from pathlib import Path

from argparse import ArgumentParser
from os.path import exists as exist
from obspy import UTCDateTime
from numpy import nan


def get_bng_average_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="%(prog)s [arguments] <Station Database>",
        description="Program to average the orientations of the seismometer " +
        "in a station database.")
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)
    parser.add_argument(
        "-v", "--verbose",
        default=2,
        type=int,
        dest="verb",
        help="Enable Level of verbose output during processing. " +
        "(0) No Output; (1) Output Event Analysis counter; " +
        "(2) Counter and results. Default 2")
    parser.add_argument(
        "--load-location",
        default="BNG_RESULTS",
        type=str,
        dest="loadloc",
        help="Specify Load destination. Default is BNG_RESULTS " +
        "(and sub-directories based on Station Name).")
    parser.add_argument(
        "--plot",
        default=False,
        action="store_true",
        dest="showplot",
        help="Plot results at end (Default False)")
    parser.add_argument(
        "--save",
        action="store_true",
        dest="saveplot",
        default=False,
        help="Set this option if you wish to save the figure. [Default " +
        "does not save figure]")
    parser.add_argument(
        "--format",
        default="png",
        dest="fmt",
        type=str,
        help="Specify format of figure. Can be any one of the valid" +
        "matplotlib formats: 'png', 'jpg', 'eps', 'pdf'. [Default 'png']")

    # Station Selection Parameters
    stparm = parser.add_argument_group(
        title="Station Selection Parameters",
        description="Parameters to select a specific station.")
    stparm.add_argument(
        "--keys",
        dest="stkeys",
        type=str,
        default="",
        help="Specify list of Station Keys in the database to process.")

    # Select QC criteria
    qcparm = parser.add_argument_group(
        title="Quality control parameters",
        description="Quality control parameters on the estimates for "+
        "calculating the average.")
    qcparm.add_argument(
        "--cc",
        dest="cc",
        type=float,
        default=0.5,
        help="Threshold for cross-correlation betwen vertical and radial "+
        "components. [Default 0.5]")
    qcparm.add_argument(
        "--snr",
        dest="snr",
        type=float,
        default=5.,
        help="Threshold for signal-to-noise ratio on vertical component, "+
        "in dB. [Default 5.]")
    qcparm.add_argument(
        "--TR",
        dest="TR",
        type=float,
        default=0.5,
        help="Threshold for transverse to radial ratio (1 - T/R). "+
        "[Default 0.5]")
    qcparm.add_argument(
        "--RZ",
        dest="RZ",
        type=float,
        default=-1.,
        help="Threshold for radial to vertical ratio (1 - R/Z). "+
        "[Default -1.]")


    # Parse Arguments
    args = parser.parse_args(argv)

    # Check inputs
    #if len(args) != 1: parser.error("Need station database file")
    # indb=args[0]
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    return args


def main(args=None):

    if args is None:
        # Run Input Parser
        args = get_bng_average_arguments()

    # Load Database
    db, stkeys = stdb.io.load_db(fname=args.indb, keys=args.stkeys)

    # Loop over station keys
    for stkey in list(stkeys):

        sta = db[stkey]

        # Input directory
        indir = Path(args.loadloc) / stkey.upper()
        if not indir.exists():
            raise(Exception("Directory does not exist: ", indir, ", aborting"))

        # Temporary print locations
        tlocs = sta.location
        if len(tlocs) == 0:
            tlocs = ['']
        for il in range(0, len(tlocs)):
            if len(tlocs[il]) == 0:
                tlocs[il] = "--"
        sta.location = tlocs

        # Update Display
        if args.verb > 1:
            print("|==============================================|")
            print("|                   {0:>8s}                   |".format(
                sta.station))
            print("|==============================================|")
            print("|  Station: {0:>2s}.{1:5s}                           |".format(
                sta.network, sta.station))
            print("|      Channel: {0:2s}; Locations: {1:15s} |".format(
                sta.channel, ",".join(tlocs)))
            print("|      Lon: {0:7.2f}; Lat: {1:6.2f}               |".format(
                sta.longitude, sta.latitude))
            print("| Input Directory: ", args.loadloc)
            print("| Plot Results: ", args.showplot)
            print("|")


        phi = []; cc = []; snr = []; TR = []; RZ = []; baz = []; mag = []
        for folder in os.listdir(indir):

            # Load meta data
            filename = indir / folder / "Meta_data.pkl"
            if not filename.is_file():
                continue
            meta = pickle.load(open(filename, 'rb'))

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

        # Set conditions for good result
        snrp = snr>args.snr
        ccp = cc>args.cc
        TRp = TR>args.TR
        RZp = RZ>args.RZ

        # Indices where conditions are met
        ind = snrp*ccp*TRp*RZp

        # Get estimate and uncertainty
        val, err = utils.estimate(phi, ind)

        # output results to termianl
        print("|    B-N-G mean, error, data included: " +
              "{0:.1f}, {1:.1f}, {2}".format(val, err, np.sum(ind)))
        print()

        if np.sum(np.isnan(np.array([val, err])))>0:
            continue

        if args.showplot or args.saveplot:

            plot = plotting.plot_bng_conditions(stkey, snr, cc, TR, RZ, ind)

            # save figure
            if args.saveplot:
                figname = indir / ('conditions.' + args.fmt)
                try:
                    plot.savefig(figname, fmt=args.fmt)
                except:
                    plot.savefig(figname, format=args.fmt)

            if args.showplot:
                plot.show()

            plot = plotting.plot_bng_results(stkey, phi, snr, cc, TR, RZ, baz, mag,
                ind, val, err)

            # save figure
            if args.saveplot:
                figname = indir / ('results.' + args.fmt)
                try:
                    plot.savefig(figname, fmt=args.fmt)
                except:
                    plot.savefig(figname, format=args.fmt)
            if args.showplot:
                plot.show()




if __name__ == "__main__":

    # Run main program
    main()
