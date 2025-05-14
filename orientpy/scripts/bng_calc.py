#!/usr/bin/env python

import stdb
import pickle
import numpy as np
import copy
# from numpy import nan

from obspy.clients.fdsn import Client as FDSN_Client
from obspy.clients.filesystem.sds import Client as SDS_Client
from obspy.taup import TauPyModel
from obspy import UTCDateTime

from orientpy import BNG, utils

from pathlib import Path
from argparse import ArgumentParser
from os.path import exists as exist


def get_bng_calc_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

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
        "-V", "--verbose",
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

    # Use local data directory
    Dtparm = parser.add_argument_group(
        title="Local Data Settings",
        description="Settings associated with defining " +
                    "and using a local data base of pre-downloaded " +
                    "day-long SAC or MSEED files.")
    Dtparm.add_argument(
        "--local-data",
        action="store",
        type=str,
        dest="localdata",
        default=None,
        help="Specify path containing " +
             "day-long sac or mseed files of data already downloaded. " +
             "If data exists for a seismogram is already present on disk, " +
             "it is selected preferentially over downloading the data " +
             "using the FDSN Client interface")
    Dtparm.add_argument(
        "--dtype",
        action="store",
        type=str,
        dest="dtype",
        default='MSEED',
        help="Specify the data archive file type, either SAC " +
             " or MSEED. Note the default behaviour is to search for " +
             "SAC files. Local archive files must have extensions of " +
             "'.SAC'  or '.MSEED'. These are case dependent, so specify " +
             "the correct case here.")

    # Server Settings
    Svparm = parser.add_argument_group(
        title="Server Settings",
        description="Settings associated with which datacenter to log into.")
    Svparm.add_argument(
        "--server-cat",
        action="store",
        type=str,
        dest="server_cat",
        default="IRIS",
        help="Catalogue server setting: Key string for recognized server " +
             "that provide `available_event_catalogs` service (one of '" +
             "AUSPASS', 'BGR', 'EARTHSCOPE', 'EIDA', 'EMSC', 'ETH', " +
             "'GEOFON', 'GEONET', 'GFZ', 'ICGC', 'IESDMC', 'INGV', 'IPGP', " +
             "'IRIS', 'IRISPH5', 'ISC', 'KNMI', 'KOERI', 'LMU', 'NCEDC', " +
             "'NIEP', 'NOA', 'NRCAN', 'ODC', 'ORFEUS', 'RASPISHAKE', " +
             "'RESIF', 'RESIFPH5', 'SCEDC', 'TEXNET', 'UIB-NORSAR', " +
             "'USGS', 'USP'). [Default 'IRIS']")
    Svparm.add_argument(
        "--server-wf",
        action="store",
        type=str,
        dest="server_wf",
        default="IRIS",
        help="Waveform server setting: Base URL of FDSN web service " +
             "compatible server (e.g. “http://service.iris.edu”) or " +
             "key string for recognized server (one of 'AUSPASS', 'BGR', " +
             "'EARTHSCOPE', 'EIDA', 'EMSC', 'ETH', 'GEOFON', 'GEONET', " +
             "'GFZ', 'ICGC', 'IESDMC', 'INGV', 'IPGP', 'IRIS', 'IRISPH5', " +
             "'ISC', 'KNMI', 'KOERI', 'LMU', 'NCEDC', 'NIEP', 'NOA', " +
             "'NRCAN', 'ODC', 'ORFEUS', 'RASPISHAKE', 'RESIF', 'RESIFPH5', " +
             "'SCEDC', 'TEXNET', 'UIB-NORSAR', 'USGS', 'USP'). " +
             "[Default 'IRIS']")
    Svparm.add_argument(
        "--user-auth",
        action="store",
        type=str,
        dest="userauth",
        default=None,
        help="Authentification Username and Password for the waveform " +
             "server (--user-auth='username:authpassword') to access " +
             "and download restricted data. [Default no user and password]")
    Svparm.add_argument(
        "--eida-token",
        action="store",
        type=str,
        dest="tokenfile",
        default=None,
        help="Token for EIDA authentication mechanism, see " +
             "http://geofon.gfz-potsdam.de/waveform/archive/auth/index.php. "
             "If a token is provided, argument --user-auth will be ignored. "
             "This mechanism is only available on select EIDA nodes. The " +
             "token can be provided in form of the PGP message as a string, " +
             "or the filename of a local file with the PGP message in it. " +
             "[Default None]")

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
        "--zcomp",
        dest="zcomp",
        type=str,
        default='Z',
        help="Specify the Vertical Component Channel Identifier. " +
             "[Default Z]")
    stparm.add_argument(
        "--coord-system",
        dest="nameconv",
        type=int,
        default=2,
        help="Coordinate system specification of instrument. " +
             "(0) Attempt Autodetect between 1 and 2; " +
             "(1) HZ, HN, HE; " +
             "(2) Left Handed: HZ, H2 90 CW H1; " +
             "(3) Right Handed: HZ, H2 90 CCW H1 " +
             "(4) Left Handed Numeric: H3, H2 90 CW H1 " +
             "[Default 2]")

    # Timing
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
        help="Enter length of time window following P arrival time " +
             "in seconds. [Default 15.]")
    Tmparm.add_argument(
        "--times",
        dest="tt",
        type=str,
        default=None,
        help="Enter window start and end times relative to predicted " +
             "P arrival time in seconds. Negative values imply start of " +
             "window before P wave arrival. [Default -2., 5.]")
    Tmparm.add_argument(
        "--sampling-rate",
        dest="new_sr",
        type=float,
        default=None,
        help="Specify new sampling rate in Hz. [Default no resampling]")

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

    # Processing Specifications
    Procparm = parser.add_argument_group(
        title="Processing Parameters",
        description="Parameters associated with BNG processing.")
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
        help="Specify corner frequencies in Hz as a list of two " +
             "floats. [Default 0.7,5.0]")
    Procparm.add_argument(
        "--plot",
        dest="showplot",
        default=False,
        action="store_true",
        help="Show processing step including raw and rotated " +
             "waveforms. [Default doesn't show plot]")

    # Parse Arguments
    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # Create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # Construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except Exception:
            parser.error(
                "Cannot construct UTCDateTime from start time: " + args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except Exception:
            parser.error(
                "Cannot construct UTCDateTime from end time: " + args.endT)
    else:
        args.endT = None

    # Parse restricted data settings
    if args.tokenfile is not None:
        args.userauth = [None, None]
    else:
        if args.userauth is not None:
            tt = args.userauth.split(':')
            if not len(tt) == 2:
                msg = ("Error: Incorrect Username and Password Strings for " +
                       "User Authentification")
                parser.error(msg)
            else:
                args.userauth = tt
        else:
            args.userauth = [None, None]

    # # Parse Local Data directories
    # if len(args.localdata) > 0:
    #     args.localdata = args.localdata.split(',')
    # else:
    #     args.localdata = []

    # # Check NoData Value
    # if args.ndval:
    #     args.ndval = 0.0
    # else:
    #     args.ndval = nan

    # Check Datatype specification
    if args.dtype.upper() not in ['MSEED', 'SAC']:
        parser.error(
            "Error: Local Data Archive must be of types 'SAC'" +
            "or MSEED. These must match the file extensions for " +
            " the archived data.")

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

    print()
    print("###########################################")
    print("#  _                              _       #")
    print("# | |__  _ __   __ _     ___ __ _| | ___  #")
    print("# | '_ \| '_ \ / _` |   / __/ _` | |/ __| #")
    print("# | |_) | | | | (_| |  | (_| (_| | | (__  #")
    print("# |_.__/|_| |_|\__, |___\___\__,_|_|\___| #")
    print("#              |___/_____|                #")
    print("#                                         #")
    print("###########################################")
    print()

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
        cat_client = FDSN_Client(
            base_url=args.server_cat)
        if args.verb > 1:
            print("      Done")

        # Establish client for waveforms
        if args.verb > 1:
            print("   Establishing Waveform Client...")

        if args.localdata is None:
            wf_client = FDSN_Client(
                base_url=args.server_wf,
                user=args.userauth[0],
                password=args.userauth[1],
                eida_token=args.tokenfile)
        else:
            wf_client = SDS_Client(
                args.localdata,
                format=args.dtype)

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
        tlocs = copy.copy(sta.location)
        if len(tlocs) == 0:
            tlocs = ['']
        for il in range(0, len(tlocs)):
            if len(tlocs[il]) == 0:
                tlocs.append("--")

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

            # get index of repeat events, save for later
            reps = np.unique(utils.catclean(cat))

        except Exception:
            raise Exception("  Fatal Error: Cannot download Catalogue")

        if args.verb > 1:
            print("|   Retrieved {0} unique events of {1}".format(
                len(cat.events)-len(reps), len(cat.events)))
            print()

        for i, ev in enumerate(cat):

            if i in reps:
                continue

            # Initialize BNGData object with station info
            bngdata = BNG(sta, zcomp=args.zcomp)

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
                    client=wf_client, new_sr=args.new_sr,
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
