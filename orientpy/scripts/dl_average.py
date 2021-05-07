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


def get_dl_average_arguments(argv=None):
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
        default="DL_RESULTS",
        type=str,
        dest="loadloc",
        help="Specify Load destination. [Default is DL_RESULTS " +
        "(and sub-directories based on Station Name)]")
    parser.add_argument(
        "--plot",
        default=False,
        action="store_true",
        dest="showplot",
        help="Plot results at end [Default False]")
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
    parser.add_argument(
        "--cc",
        default=0.8,
        type=float,
        dest="cc",
        help="Cross-correlation threshold for final estimate. [Default 0.8]")

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
        args = get_dl_average_arguments()

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


        R1phi = []; R1cc = []; R2phi = []; R2cc = []
        for folder in os.listdir(indir):

            # Load meta data
            filename = indir / folder / "Meta_data.pkl"
            if not filename.is_file():
                continue
            meta = pickle.load(open(filename, 'rb'))

            R1phi.append(meta.R1phi)
            R2phi.append(meta.R2phi)
            R1cc.append(meta.R1cc)
            R2cc.append(meta.R2cc)

        R1phi = np.array(R1phi).flatten()
        R1cc = np.array(R1cc).flatten()
        R2phi = np.array(R2phi).flatten()
        R2cc = np.array(R2cc).flatten()

        phi = np.concatenate((R1phi, R2phi), axis=None)
        cc = np.concatenate((R1cc, R2cc), axis=None)
        ind = cc > args.cc

        val, err = utils.estimate(phi, ind)

        # output results to termianl
        print("|    D-L mean, error, data included: " +
              "{0:.2f}, {1:.2f}, {2}".format(val, err, np.sum(ind)))
        print("|    D-L CC level: {0:.1f}".format(args.cc))
        print()

        if np.sum(np.isnan(np.array([val, err])))>0:
            continue

        if args.showplot or args.saveplot:

            plot = plotting.plot_dl_results(stkey, R1phi, R1cc, R2phi, R2cc, ind,
                val, err, phi, cc, args.cc)

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
