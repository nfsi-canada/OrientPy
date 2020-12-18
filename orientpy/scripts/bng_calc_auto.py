#!/usr/bin/env python

import stdb
import pickle
import os.path
import numpy as np
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from orientpy import BNG
from pathlib import Path

from argparse import ArgumentParser
from os.path import exists as exist
from obspy import UTCDateTime
from numpy import nan


def get_bng_calc_arguments(argv=None):
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
        default="BNG_RESULTS",
        type=str,
        dest="saveloc",
        help="Specify Save destination. Default is BNG_RESULTS " +
        "(and sub-directories based on Station Name).")
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
        dest="wlen",
        type=float,
        default=15.,
        help="Enter length of time window following P arrival time in "+
        "seconds. [Default 15.]")
    Tmparm.add_argument(
        "--times",
        dest="tt",
        type=str,
        default=None,
        help="Enter window start and end times relative to predicted P "+
        "arrival time in seconds. Negative values imply start of window "+
        "before P wave arrival. [Default -2., 5.]")

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
        "--max-mag",
        dest="maxmag",
        type=float,
        default=9.,
        help="Specify the maximum magnitude of Earthquakes to use in " +
        "the catalogue search. [Default 9.]")
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
        default=1000.,
        help="Specify maximum Earthquake Depth (km). [Default no limit]")
    Eqparm.add_argument(
        "--discard-catalogue",
        dest="savecat",
        default=True,
        action="store_false",
        help="Specify to discard the eq catalogue after processing.")

    # Processing Specifications
    Procparm = parser.add_argument_group(
        title="Processing Parameters",
        description="Parameters associated with BNG processing.")
    Procparm.add_argument(
        "--new-sampling-rate",
        dest="new_sr",
        type=float,
        default=None,
        help="Specify new sampling rate in Hz. [Default no resampling]")
    Procparm.add_argument(
        "--dphi",
        dest="dphi",
        type=float,
        default=0.1,
        help="Specify angle interval for search, in degrees. [Default 0.1]")
    Procparm.add_argument(
        "--bp",
        dest="bp",
        type=str,
        default=None,
        help="Specify corner frequencies in Hz as a list of two floats. "+
        "[Default 0.7,5.0]")
    Procparm.add_argument(
        "--plot",
        dest="showplot",
        default=False,
        action="store_true",
        help="Show processing step including raw and rotated waveforms. "+
        "[Default doesn't show plot]")

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

    if args.bp is not None:
        args.bp = [float(val) for val in args.bp.split(',')]
        args.bp = sorted(args.bp)
        if (len(args.bp)) != 2:
            parser.error(
                "Error: --bp should contain 2 " +
                "comma-separated floats")

    if args.tt is not None:
        args.tt = [float(val) for val in args.tt.split(',')]
        args.tt = sorted(args.tt)
        if (len(args.tt)) != 2:
            parser.error(
                "Error: --times should contain 2 " +
                "comma-separated floats")
    else:
        args.tt = [-2., 5.]

    return args


def main(args=None):

    if args is None:
        # Run Input Parser
        args = get_bng_calc_arguments()

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
            print("|   Mag:   >{0:3.1f}".format(args.minmag) +
                  "                                |")
            print("|   Min Distance:  {0:.1f}".format(args.mindist))
            print("|   Max Distance:  {0:.1f}".format(args.maxdist))
            print("|   Max Depth:     {0:.1f}".format(args.maxdep))

        # Retrieve Event Catalogue
        if args.verb > 1:
            print("|   Request Event Catalogue...                 |")
            print("| ...                                          |")

        try:
            cat = cat_client.get_events(
                starttime=tstart, endtime=tend,
                minmagnitude=args.minmag, maxmagnitude=args.maxmag)
        except:
            raise(Exception("  Fatal Error: Cannot download Catalogue"))

        if args.verb > 1:
            print("|   Retrieved {0} events ".format(len(cat.events)))
            print()

        for i, ev in enumerate(cat):

            # Initialize BNGData object with station info
            bngdata = BNG(sta)

            # Add event to object
            accept = bngdata.add_event(
                ev, gacmin=args.mindist, gacmax=args.maxdist,
                depmax=args.maxdep, returned=True)

            # Define time stamp
            yr = str(bngdata.meta.time.year).zfill(4)
            jd = str(bngdata.meta.time.julday).zfill(3)
            hr = str(bngdata.meta.time.hour).zfill(2)

            # If event is accepted (data exists)
            if accept:

                # Get travel time info
                tpmodel = TauPyModel(model='iasp91')

                # Get Travel times
                arrivals = tpmodel.get_travel_times(
                    distance_in_degree=bngdata.meta.gac,
                    source_depth_in_km=bngdata.meta.dep,
                    phase_list=['P', 'PP'])

                # Get first P wave arrival among P and PP
                arrival = arrivals[0]

                # Attributes from parameters
                bngdata.meta.ttime = arrival.time
                bngdata.meta.phase = arrival.name

                # Display Event Info
                print(" ")
                print("**************************************************")
                print("* ({0:d}/{1:d}):  {2:13s} {3}".format(
                    i+1, len(cat), bngdata.meta.time.strftime(
                        "%Y%m%d_%H%M%S"), stkey))
                if args.verb > 1:
                    print("*   Phase: {}".format(bngdata.meta.phase))
                    print("*   Origin Time: " +
                          bngdata.meta.time.strftime("%Y-%m-%d %H:%M:%S"))
                    print(
                        "*   Lat: {0:6.2f};        Lon: {1:7.2f}".format(
                            bngdata.meta.lat, bngdata.meta.lon))
                    print(
                        "*   Dep: {0:6.2f} km;     Mag: {1:3.1f}".format(
                            bngdata.meta.dep, bngdata.meta.mag))
                    print("*   Dist: {0:7.2f} km;".format(bngdata.meta.epi_dist) +
                          "   Epi dist: {0:6.2f} deg\n".format(bngdata.meta.gac) +
                          "*   Baz:  {0:6.2f} deg;".format(bngdata.meta.baz) +
                          "   Az: {0:6.2f} deg".format(bngdata.meta.az))

                # Event Folder
                timekey = bngdata.meta.time.strftime("%Y%m%d_%H%M%S")
                evtdir = outdir / timekey
                evtdata = evtdir / 'Raw_data.pkl'
                evtmeta = evtdir / 'Meta_data.pkl'

                # Check if BNG data already exist and overwrite has been set
                if evtdir.exists():
                    if evtdata.exists():
                        if not args.ovr:
                            continue

                # Get data
                t1 = arrival.time - args.wlen
                t2 = arrival.time + args.wlen
                has_data = bngdata.download_data(
                    client=wf_client, stdata=args.localdata,
                    ndval=args.ndval, new_sr=args.new_sr,
                    t1=t1, t2=t2, returned=True, verbose=args.verb)

                if not has_data:
                    continue

                # Create Folder if it doesn't exist
                if not evtdir.exists():
                    evtdir.mkdir(parents=True)

                # Save raw Traces
                pickle.dump(bngdata.data, open(evtdata, "wb"))

                # Calculate BNG orientation
                bngdata.calc(args.dphi, args.wlen, args.tt,
                    bp=args.bp, showplot=args.showplot)

                if args.verb > 1:
                    print("* PHI: {}".format(bngdata.meta.phi))
                    print("* SNR: {}".format(bngdata.meta.snr))
                    print("* CC: {}".format(bngdata.meta.cc))
                    print("* 1-T/R: {}".format(bngdata.meta.TR))
                    print("* 1-R/Z: {}".format(bngdata.meta.RZ))

                # Save event meta data
                pickle.dump(bngdata.meta, open(evtmeta, "wb"))


if __name__ == "__main__":

    # Run main program
    main()
