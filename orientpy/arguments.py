'''
Utilities for orient module

AJS September 2017
'''

from argparse import ArgumentParser
from os.path import exists as exist
from obspy import UTCDateTime
from numpy import nan


def get_bng_calc_arguments():
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <Station Database>",
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
        help="Specify Save destination. Default is DLOPY_RESULTS " +
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
    args = parser.parse_args()

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
            parser.errer(
                "Error: Incorrect Username and Password Strings " +
                "for User Authentification")
        else:
            args.UserAuth = tt
    else:
        args.UserAuth = []

    # #-- Check existing file behaviour
    # if opts.skip and opts.ovr:
    # opts.skip=False
    # opts.ovr=False

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

    if args.tt is None:
        args.tt = [-2., 5.]
    else:
        args.tt = [float(val) for val in args.tt.split(',')]
        args.tt = sorted(args.tt)
        if (len(args.tt)) != 2:
            parser.error(
                "Error: --times should contain 2 " +
                "comma-separated floats")

    return args


def get_bng_average_arguments():
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <Station Database>",
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
        help="Specify Load destination. Default is DLOPY_RESULTS " +
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
    args = parser.parse_args()

    # Check inputs
    #if len(args) != 1: parser.error("Need station database file")
    # indb=args[0]
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    return args


def get_dl_calc_arguments():
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <Station Database>",
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
        action="store_false",
        dest="ovr",
        help="Overwrite existing data on disk. [Default False]")
    parser.add_argument(
        "--save-location",
        default="DL_RESULTS",
        type=str,
        dest="saveloc",
        help="Specify Save destination. Default is DLOPY_RESULTS " +
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
    args = parser.parse_args()

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
            parser.errer(
                "Error: Incorrect Username and Password Strings " +
                "for User Authentification")
        else:
            args.UserAuth = tt
    else:
        args.UserAuth = []

    # #-- Check existing file behaviour
    # if opts.skip and opts.ovr:
    # opts.skip=False
    # opts.ovr=False

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


def get_dl_average_arguments():
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <Station Database>",
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
        help="Specify Load destination. Default is DLOPY_RESULTS " +
        "(and sub-directories based on Station Name).")
    parser.add_argument(
        "--plot",
        default=False,
        action="store_true",
        dest="showplot",
        help="Plot results at end (Default False)")
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
    args = parser.parse_args()

    # Check inputs
    #if len(args) != 1: parser.error("Need station database file")
    # indb=args[0]
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')


    return args

