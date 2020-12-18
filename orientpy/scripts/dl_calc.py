#!/usr/bin/env python

# script to test the orient.py running script.

import stdb
import pickle
import numpy as np
from obspy.clients.fdsn import Client
from orientpy import DL, utils
from pathlib import Path

from argparse import ArgumentParser
from os.path import exists as exist
from obspy import UTCDateTime
from numpy import nan


def get_dl_calc_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="%(prog)s [arguments] <Station Database>",
        description="Program to compute the orientation of the components " +
        "of a station based on those in a station database.")
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
        "-O", "--overwrite",
        default=False,
        action="store_true",
        dest="ovr",
        help="Overwrite existing data on disk. [Default False]")
    parser.add_argument(
        "--save-location",
        default="DL_RESULTS",
        type=str,
        dest="saveloc",
        help="Specify Save destination. [Default is DL_RESULTS " +
        "(and sub-directories based on Station Name)]")
    parser.add_argument(
        "--no-save-progress",
        default=True,
        action="store_false",
        dest="constsave",
        help="Do not save progress during processing.")

    # Use local data directory
    Dtparm = parser.add_argument_group(
        title="Local Data Settings",
        description="Settings associated with defining and using a " +
        "local data base of pre-downloaded day-long SAC files.")
    Dtparm.add_argument(
        "--local-data",
        action="store",
        type=str,
        dest="localdata",
        default="",
        help="Specify a comma separated list of paths containing " +
        "day-long sac files of data already downloaded. If data exists " +
        "for a seismogram is already present on disk, it is selected " +
        "preferentially over downloading the data using the Client interface")
    Dtparm.add_argument(
        "--no-data-zero",
        action="store_true",
        dest="ndval",
        default=False,
        help="Specify to force missing data to be set as zero, rather " +
        "than default behaviour. [Default sets to nan]")
    Dtparm.add_argument(
        "--no-local-net",
        action="store_false",
        dest="useNet",
        default=True,
        help="Specify to prevent using the Network code in the search " +
        "for local data (sometimes for CN stations the dictionary name " +
        "for a station may disagree with that in the filename. " +
        "[Default Network used]")

    # Server Settings
    Svparm = parser.add_argument_group(
        title="Server Settings",
        description="Settings associated with which datacenter to log into.")
    Svparm.add_argument(
        "--catalogue-source",
        action="store",
        type=str,
        dest="cat_client",
        default="IRIS",
        help="Specify the server to connect to for the event catalogue. " +
        "Options include: BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, " +
        "LMU, NCEDC, NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. " +
        "[Default IRIS]")
    Svparm.add_argument(
        "--waveform-source",
        action="store",
        type=str,
        dest="wf_client",
        default="IRIS",
        help="Specify the server to connect to for the waveform data. " +
        "Options include: BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, " +
        "LMU, NCEDC, NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. " +
        "[Default IRIS]")
    Svparm.add_argument(
        "-U",
        "--User-Auth",
        action="store",
        type=str,
        dest="UserAuth",
        default="",
        help="Enter your Authentification Username and Password for the " +
        "waveform server (--User-Auth='username:authpassword') to access " +
        "and download restricted data. [Default no user and password]")

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
    stparm.add_argument(
        "-c", "--coord-system",
        dest="nameconv",
        type=int,
        default=2,
        help="Coordinate system specification of instrument. " +
        "(0) Attempt Autodetect between 1 and 2; (1) HZ, HN, HE; " +
        "(2) Left Handed: HZ, H2 90 CW H1; (3) Right Handed: HZ, H2 90 CCW " +
        "H1. [Default 2]")

    #-- Timing
    Tmparm = parser.add_argument_group(
        title="Timing Parameters",
        description="Parameters associated with event timing and window " +
        "length.")
    Tmparm.add_argument(
        "--start",
        dest="startT",
        type=str,
        default="",
        help="Enter Start date for event catalogue search. Note, more " +
        "recent of this value or station start date will be used.")
    Tmparm.add_argument(
        "--end",
        dest="endT",
        type=str,
        default="",
        help="Enter End date for event catalogue search. Note, less " +
        "recent of this or the station end date will be used.")
    Tmparm.add_argument(
        "--window",
        dest="twin",
        type=int,
        default=0,
        help="Enter time window length in days. A non-zero value will " +
        "cause the results to repeat for each set of twin days in " +
        "the operating window, calculating the change in orientation over " +
        "time. [Default 0]")

    # EQ Specifications
    Eqparm = parser.add_argument_group(
        title="Earthquake Selection Criteria",
        description="Parameters associated with selecing the subset of " +
        "earthquakes to use in calculations.")
    Eqparm.add_argument(
        "--min-mag",
        dest="minmag",
        type=float,
        default=5.5,
        help="Specify the minimum magnitude of Earthquakes to use in " +
        "the catalogue search. [Default 5.5]")
    Eqparm.add_argument(
        "--min-dist",
        dest="mindist",
        type=float,
        default=5.,
        help="Specify the minimum earthquake distance (in degrees). " +
        "[Default 5.]")
    Eqparm.add_argument(
        "--max-dist",
        dest="maxdist",
        type=float,
        default=175.,
        help="Specify the maximum earthquake distance (in degrees). " +
        "[Default 175.]")
    Eqparm.add_argument(
        "--max-dep",
        dest="maxdep",
        type=float,
        default=150.,
        help="Specify maximum Earthquake Depth (km). [Default 150.]")
    Eqparm.add_argument(
        "--discard-catalogue",
        dest="savecat",
        default=True,
        action="store_false",
        help="Specify to discard the eq catalogue after processing.")

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

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Cannot construct UTCDateTime from start time: " + args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Cannot construct UTCDateTime from end time: " + args.endT)
    else:
        args.endT = None

    # Parse User Authentification
    if not len(args.UserAuth) == 0:
        tt = args.UserAuth.split(':')
        if not len(tt) == 2:
            parser.error(
                "Error: Incorrect Username and Password Strings " +
                "for User Authentification")
        else:
            args.UserAuth = tt
    else:
        args.UserAuth = []

    # Parse Local Data directories
    if len(args.localdata) > 0:
        args.localdata = args.localdata.split(',')
    else:
        args.localdata = []

    # Check NoData Value
    if args.ndval:
        args.ndval = 0.0
    else:
        args.ndval = nan

    return args


def main(args=None):

    if args is None:
        # Run Input Parser
        args = get_dl_calc_arguments()

    # Load Database
    db, stkeys = stdb.io.load_db(fname=args.indb, keys=args.stkeys)

    # Loop over station keys
    for stkey in list(stkeys):

        sta = db[stkey]

        # Output directory
        outdir = Path(args.saveloc) / Path(stkey.upper())
        if not outdir.exists():
            outdir.mkdir(parents=True)

        # Establish client for catalogue
        if args.verb > 1:
            print("   Establishing Catalogue Client...")
        cat_client = Client(args.cat_client)
        if args.verb > 1:
            print("      Done")

        # Establish client for waveforms
        if args.verb > 1:
            print("   Establishing Waveform Client...")
        if len(args.UserAuth) == 0:
            wf_client = Client(args.wf_client)
        else:
            wf_client = Client(args.wf_client,
                               user=args.UserAuth[0],
                               password=args.UserAuth[1])
        if args.verb > 1:
            print("      Done")
            print(" ")

        # Get catalogue search start time
        if args.startT is None:
            tstart = sta.startdate
        else:
            tstart = args.startT

        # Get catalogue search end time
        if args.endT is None:
            tend = sta.enddate
        else:
            tend = args.endT
        if tstart > sta.enddate or tend < sta.startdate:
            continue

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
            print("|      Start time: {0:19s}         |".format(
                sta.startdate.strftime("%Y-%m-%d %H:%M:%S")))
            print("|      End time:   {0:19s}         |".format(
                sta.enddate.strftime("%Y-%m-%d %H:%M:%S")))
            print("| Output Directory: ", args.saveloc)
            print("| Save Progress: ", args.constsave)
            print("|----------------------------------------------|")
            print("| Searching Possible events:                   |")
            print("|   Start: {0:19s}                 |".format(
                tstart.strftime("%Y-%m-%d %H:%M:%S")))
            print("|   End:   {0:19s}                 |".format(
                tend.strftime("%Y-%m-%d %H:%M:%S")))
            print("|   Mag:   >{0:3.1f}", format(args.minmag) +
                  "                       |")
            print("|   Min Distance:  {0:.1f}".format(args.mindist))
            print("|   Max Distance:  {0:.1f}".format(args.maxdist))
            print("|   Max Depth:     {0:.1f}".format(args.maxdep))

        # Retrieve Event Catalogue
        if args.verb > 1:
            print("|   Request Event Catalogue...                 |")
            print("| ...                                          |")

        try:
            cat = cat_client.get_events(starttime=tstart, endtime=tend,
                                        minmagnitude=args.minmag)

            # get index of repeat events, save for later
            reps = np.unique(utils.catclean(cat))

        except:
            raise(Exception("  Fatal Error: Cannot download Catalogue"))

        if args.verb > 1:
            print("|   Retrieved {0} events ".format(len(cat.events)))
            print()

        for i, ev in enumerate(cat):

            if i in reps:
                continue

            # Initialize BNGData object with station info
            dldata = DL(sta)

            # Add event to object
            accept = dldata.add_event(
                ev, gacmin=args.mindist, gacmax=args.maxdist,
                depmax=args.maxdep, returned=True)

            # Define time stamp
            yr = str(dldata.meta.time.year).zfill(4)
            jd = str(dldata.meta.time.julday).zfill(3)
            hr = str(dldata.meta.time.hour).zfill(2)

            # If event is accepted (data exists)
            if accept:

                # Display Event Info
                print(" ")
                print("**************************************************")
                print("* ({0:d}/{1:d}):  {2:13s} {3}".format(
                    i+1, len(cat), dldata.meta.time.strftime(
                        "%Y%m%d_%H%M%S"), stkey))
                if args.verb > 1:
                    print("*   Origin Time: " +
                          dldata.meta.time.strftime("%Y-%m-%d %H:%M:%S"))
                    print(
                        "*   Lat: {0:6.2f};        Lon: {1:7.2f}".format(
                            dldata.meta.lat, dldata.meta.lon))
                    print(
                        "*   Dep: {0:6.2f} km;     Mag: {1:3.1f}".format(
                            dldata.meta.dep, dldata.meta.mag))
                    print("*   Dist: {0:7.2f} km;".format(dldata.meta.epi_dist) +
                          "   Epi dist: {0:6.2f} deg\n".format(dldata.meta.gac) +
                          "*   Baz:  {0:6.2f} deg;".format(dldata.meta.baz) +
                          "   Az: {0:6.2f} deg".format(dldata.meta.az))

                # Event Folder
                timekey = dldata.meta.time.strftime("%Y%m%d_%H%M%S")
                evtdir = outdir / timekey
                evtdata = evtdir / 'Raw_data.pkl'
                evtmeta = evtdir / 'Meta_data.pkl'

                # Check if DL data already exist and overwrite has been set
                if evtdir.exists():
                    if evtdata.exists():
                        if not args.ovr:
                            continue

                # Get data
                t1 = 0.
                t2 = 4.*60.*60.
                has_data = dldata.download_data(
                    client=wf_client, stdata=args.localdata,
                    ndval=args.ndval, new_sr=2., t1=t1, t2=t2,
                    returned=True, verbose=args.verb)

                if not has_data:
                    continue

                # Check data length
                if utils.checklen(dldata.data, 4.*60.*60.):
                    print("      Error: Length Incorrect")
                    continue

                # Create Folder if it doesn't exist
                if not evtdir.exists():
                    evtdir.mkdir(parents=True)

                # Save raw Traces
                pickle.dump(dldata.data, open(evtdata, "wb"))

                # Calculate DL orientation
                dldata.calc(showplot=False)

                if args.verb > 1:
                    print("* R1PHI: {}".format(dldata.meta.R1phi))
                    print("* R2PHI: {}".format(dldata.meta.R2phi))
                    print("* R1CC: {}".format(dldata.meta.R1cc))
                    print("* R2CC: {}".format(dldata.meta.R2cc))

                # Save event meta data
                pickle.dump(dldata.meta, open(evtmeta, "wb"))


if __name__ == "__main__":

    # Run main program
    main()
