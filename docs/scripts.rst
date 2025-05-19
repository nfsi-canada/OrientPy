Scripts
=======

The program :mod:`~orientpy` is meant to be run using command-line scripts that
are described below. These can be used in bash scripts to automate data processing. 
These scripts use classes defined in :mod:`~orientpy` to process single-station
and single-event seismograms, which are then aggregated to generate a single 
estimate of station orientation per method. There are three methods with accompanying
scripts that can be used to determine station orientation, which are described below. 
All of them use a station database provided as a :class:`~stdb.StDb` dictionary. 


BNG
+++

``bng_calc``
*****************

Description
-----------

Downloads three-component seismograms based on a catalogue of earthquakes 
and performs automated estimation of P-wave polarization. Station selection is 
specified by a network and station code. The database is provided as a 
:class:`~stdb.StDb` dictionary. This method can be used with teleseismic P-wave
data or regional earthquake data, by specifying the appropriate options accordingly.
Each usable station-event seismograms are used to calculate an estimate of station
orientation. For each estimate, a number of quality-control parameters are calculated
to help in the aggregation to produce a final estimate of station orientation.

Usage
-----

.. code-block::

    $ bng_calc -h

    ###########################################
    #  _                              _       #
    # | |__  _ __   __ _     ___ __ _| | ___  #
    # | '_ \| '_ \ / _` |   / __/ _` | |/ __| #
    # | |_) | | | | (_| |  | (_| (_| | | (__  #
    # |_.__/|_| |_|\__, |___\___\__,_|_|\___| #
    #              |___/_____|                #
    #                                         #
    ###########################################

    usage: bng_calc [arguments] <Station Database>

    Program to compute the orientation of the components of a station based on
    those in a station database.

    positional arguments:
      indb                  Station Database to process from.

    optional arguments:
      -h, --help            show this help message and exit
      -V VERB, --verbose VERB
                            Enable Level of verbose output during processing. (0)
                            No Output; (1) Output Event Analysis counter; (2)
                            Counter and results. Default 2
      -O, --overwrite       Overwrite existing data on disk. [Default False]
      --save-location SAVELOC
                            Specify Save destination. Default is BNG_RESULTS (and
                            sub-directories based on Station Name).

    Local Data Settings:
      Settings associated with defining and using a local data base of pre-
      downloaded day-long SAC files.

      --local-data LOCALDATA
                            Specify absolute path to a SeisComP Data Structure (SDS) archive
                            containing day-long SAC or MSEED files(e.g., --local-
                            data=/Home/username/Data/SDS). See
                            https://www.seiscomp.de/seiscomp3/doc/applications/slarchive/SDS.html
                            for details on the SDS format. If this option is used, it takes
                            precedence over the --server-wf settings.
      --dtype DTYPE         Specify the data archive file type, either SAC or MSEED. Note the
                            default behaviour is to search for SAC files. Local archive files
                            must have extensions of '.SAC' or '.MSEED'. These are case dependent,
                            so specify the correct case here.

    Server Settings:
      Settings associated with which datacenter to log into.

      --server-cat SERVER_CAT
                            Catalogue server setting: Key string for recognized
                            server that provide `available_event_catalogs` service
                            (one of 'AUSPASS', 'BGR', 'EARTHSCOPE', 'EIDA',
                            'EMSC', 'ETH', 'GEOFON', 'GEONET', 'GFZ', 'ICGC',
                            'IESDMC', 'INGV', 'IPGP', 'IRIS', 'IRISPH5', 'ISC',
                            'KNMI', 'KOERI', 'LMU', 'NCEDC', 'NIEP', 'NOA',
                            'NRCAN', 'ODC', 'ORFEUS', 'RASPISHAKE', 'RESIF',
                            'RESIFPH5', 'SCEDC', 'TEXNET', 'UIB-NORSAR', 'USGS',
                            'USP'). [Default 'IRIS']
      --server-wf SERVER_WF
                            Waveform server setting: Base URL of FDSN web service
                            compatible server (e.g. “http://service.iris.edu”) or
                            key string for recognized server (one of 'AUSPASS',
                            'BGR', 'EARTHSCOPE', 'EIDA', 'EMSC', 'ETH', 'GEOFON',
                            'GEONET', 'GFZ', 'ICGC', 'IESDMC', 'INGV', 'IPGP',
                            'IRIS', 'IRISPH5', 'ISC', 'KNMI', 'KOERI', 'LMU',
                            'NCEDC', 'NIEP', 'NOA', 'NRCAN', 'ODC', 'ORFEUS',
                            'RASPISHAKE', 'RESIF', 'RESIFPH5', 'SCEDC', 'TEXNET',
                            'UIB-NORSAR', 'USGS', 'USP'). [Default 'IRIS']
      --user-auth USERAUTH  Authentification Username and Password for the
                            waveform server (--user-auth='username:authpassword')
                            to access and download restricted data. [Default no
                            user and password]

      --eida-token TOKENFILE
                            Token for EIDA authentication mechanism, see
                            http://geofon.gfz-
                            potsdam.de/waveform/archive/auth/index.php. If a token
                            is provided, argument --user-auth will be ignored.
                            This mechanism is only available on select EIDA nodes.
                            The token can be provided in form of the PGP message
                            as a string, or the filename of a local file with the
                            PGP message in it. [Default None]

    Station Selection Parameters:
      Parameters to select a specific station.

      --keys STKEYS         Specify list of Station Keys in the database to
                            process.
      --zcomp ZCOMP         Specify the Vertical Component Channel Identifier.
                            [Default Z].
      -c NAMECONV, --coord-system NAMECONV
                            Coordinate system specification of instrument. (0)
                            Attempt Autodetect between 1 and 2; (1) HZ, HN, HE;
                            (2) Left Handed: HZ, H2 90 CW H1; (3) Right Handed:
                            HZ, H2 90 CCW H1. **Note**: this option
                            is not yet implemented. [Default 2]

    Timing Parameters:
      Parameters associated with event timing and window length.

      --start STARTT        Enter Start date for event catalogue search. Note,
                            more recent of this value or station start date will
                            be used.
      --end ENDT            Enter End date for event catalogue search. Note, less
                            recent of this or the station end date will be used.
      --window WLEN         Enter length of time window following P arrival time
                            in seconds. [Default 15.]
      --times TT            Enter window start and end times relative to predicted
                            P arrival time in seconds. Negative values imply start
                            of window before P wave arrival. [Default -2., 5.]

    Earthquake Selection Criteria:
      Parameters associated with selecing the subset of earthquakes to use in
      calculations.

      --min-mag MINMAG      Specify the minimum magnitude of Earthquakes to use in
                            the catalogue search. [Default 5.5]
      --max-mag MAXMAG      Specify the maximum magnitude of Earthquakes to use in
                            the catalogue search. [Default 9.]
      --min-dist MINDIST    Specify the minimum earthquake distance (in degrees).
                            [Default 5.]
      --max-dist MAXDIST    Specify the maximum earthquake distance (in degrees).
                            [Default 175.]
      --max-dep MAXDEP      Specify maximum Earthquake Depth (km). [Default no
                            limit]

    Processing Parameters:
      Parameters associated with BNG processing.

      --new-sampling-rate NEW_SR
                            Specify new sampling rate in Hz. [Default no
                            resampling]
      --dphi DPHI           Specify angle interval for search, in degrees.
                            [Default 0.1]
      --bp BP               Specify corner frequencies in Hz as a list of two
                            floats. [Default 0.7,5.0]
      --plot                Show processing step including raw and rotated
                            waveforms. [Default doesn't show plot]


``bng_average``
***************

Description
-----------

Collects the estimated azimuths previously calculated and calculates the
mean value after some quality control thresholding based on the rotated 
waveforms. The error is obtained from a bootstrap analysis of robust estimates.

Usage
-----

.. code-block::

    $ bng_average -h

    ###############################################################
    #  _                                                          #
    # | |__  _ __   __ _     __ ___   _____ _ __ __ _  __ _  ___  #
    # | '_ \| '_ \ / _` |   / _` \ \ / / _ \ '__/ _` |/ _` |/ _ \ #
    # | |_) | | | | (_| |  | (_| |\ V /  __/ | | (_| | (_| |  __/ #
    # |_.__/|_| |_|\__, |___\__,_| \_/ \___|_|  \__,_|\__, |\___| #
    #              |___/_____|                        |___/       #
    #                                                             #
    ###############################################################

    usage: bng_average [arguments] <Station Database>

    Program to average the orientations of the seismometer in a station database.

    positional arguments:
      indb                  Station Database to process from.

    optional arguments:
      -h, --help            show this help message and exit
      -V VERB, --verbose VERB
                            Enable Level of verbose output during processing. (0)
                            No Output; (1) Output Event Analysis counter; (2)
                            Counter and results. Default 2
      --load-location LOADLOC
                            Specify Load destination. Default is BNG_RESULTS (and
                            sub-directories based on Station Name).
      --plot                Plot results at end (Default False)
      --save                Set this option if you wish to save the figure.
                            [Default does not save figure]
      --format FMT          Specify format of figure. Can be any one of the
                            validmatplotlib formats: 'png', 'jpg', 'eps', 'pdf'.
                            [Default 'png']

    Station Selection Parameters:
      Parameters to select a specific station.

      --keys STKEYS         Specify list of Station Keys in the database to
                            process.

    Quality control parameters:
      Quality control parameters on the estimates for calculating the average.

      --cc CC               Threshold for cross-correlation betwen vertical and
                            radial components. [Default 0.5]
      --snr SNR             Threshold for signal-to-noise ratio on vertical
                            component, in dB. [Default 5.]
      --TR TR               Threshold for transverse to radial ratio (1 - T/R).
                            [Default 0.5]
      --RZ RZ               Threshold for radial to vertical ratio (1 - R/Z).
                            [Default -1.]

DL
++

``dl_calc``
***********

Description
-----------

Downloads three-component seismograms based on a catalogue of earthquakes 
and performs automated estimation of Rayleigh-wave polarization at a number of
periods and for the direct and complementary globe-encircling path. Station 
selection is specified by a network and station code. The database is provided as 
a :class:`~stdb.StDb` dictionary. Each usable station-event seismograms are used 
to calculate an estimate of station orientation. For each estimate, the 
cross-correlation between the radial and Hilbert-transformed vertical components
is calculated and is used later in selecting which estimates are used in the final
estimate of station orientation.

Usage
-----

.. code-block::

    $ dl_calc -h

    #################################
    #      _ _              _       #
    #   __| | |    ___ __ _| | ___  #
    #  / _` | |   / __/ _` | |/ __| #
    # | (_| | |  | (_| (_| | | (__  #
    #  \__,_|_|___\___\__,_|_|\___| #
    #        |_____|                #
    #                               #
    #################################

    usage: dl_calc [arguments] <Station Database>

    Program to compute the orientation of the components of a station based on
    those in a station database.

    positional arguments:
      indb                  Station Database to process from.

    optional arguments:
      -h, --help            show this help message and exit
      -V VERB, --verbose VERB
                            Enable Level of verbose output during processing. (0)
                            No Output; (1) Output Event Analysis counter; (2)
                            Counter and results. Default 2
      -O, --overwrite       Overwrite existing data on disk. [Default False]
      --save-location SAVELOC
                            Specify Save destination. [Default is DL_RESULTS (and
                            sub-directories based on Station Name)]

    Local Data Settings:
      Settings associated with defining and using a local data base of pre-
      downloaded day-long SAC files.

      --local-data LOCALDATA
                            Specify absolute path to a SeisComP Data Structure (SDS) archive
                            containing day-long SAC or MSEED files(e.g., --local-
                            data=/Home/username/Data/SDS). See
                            https://www.seiscomp.de/seiscomp3/doc/applications/slarchive/SDS.html
                            for details on the SDS format. If this option is used, it takes
                            precedence over the --server-wf settings.
      --dtype DTYPE         Specify the data archive file type, either SAC or MSEED. Note the
                            default behaviour is to search for SAC files. Local archive files
                            must have extensions of '.SAC' or '.MSEED'. These are case dependent,
                            so specify the correct case here.

    Server Settings:
      Settings associated with which datacenter to log into.

      --server-cat SERVER_CAT
                            Catalogue server setting: Key string for recognized
                            server that provide `available_event_catalogs` service
                            (one of 'AUSPASS', 'BGR', 'EARTHSCOPE', 'EIDA',
                            'EMSC', 'ETH', 'GEOFON', 'GEONET', 'GFZ', 'ICGC',
                            'IESDMC', 'INGV', 'IPGP', 'IRIS', 'IRISPH5', 'ISC',
                            'KNMI', 'KOERI', 'LMU', 'NCEDC', 'NIEP', 'NOA',
                            'NRCAN', 'ODC', 'ORFEUS', 'RASPISHAKE', 'RESIF',
                            'RESIFPH5', 'SCEDC', 'TEXNET', 'UIB-NORSAR', 'USGS',
                            'USP'). [Default 'IRIS']
      --server-wf SERVER_WF
                            Waveform server setting: Base URL of FDSN web service
                            compatible server (e.g. “http://service.iris.edu”) or
                            key string for recognized server (one of 'AUSPASS',
                            'BGR', 'EARTHSCOPE', 'EIDA', 'EMSC', 'ETH', 'GEOFON',
                            'GEONET', 'GFZ', 'ICGC', 'IESDMC', 'INGV', 'IPGP',
                            'IRIS', 'IRISPH5', 'ISC', 'KNMI', 'KOERI', 'LMU',
                            'NCEDC', 'NIEP', 'NOA', 'NRCAN', 'ODC', 'ORFEUS',
                            'RASPISHAKE', 'RESIF', 'RESIFPH5', 'SCEDC', 'TEXNET',
                            'UIB-NORSAR', 'USGS', 'USP'). [Default 'IRIS']
      --user-auth USERAUTH  Authentification Username and Password for the
                            waveform server (--user-auth='username:authpassword')
                            to access and download restricted data. [Default no
                            user and password]
      --eida-token TOKENFILE
                            Token for EIDA authentication mechanism, see
                            http://geofon.gfz-
                            potsdam.de/waveform/archive/auth/index.php. If a token
                            is provided, argument --user-auth will be ignored.
                            This mechanism is only available on select EIDA nodes.
                            The token can be provided in form of the PGP message
                            as a string, or the filename of a local file with the
                            PGP message in it. [Default None]

    Station Selection Parameters:
      Parameters to select a specific station.

      --keys STKEYS         Specify list of Station Keys in the database to
                            process.
      --zcomp ZCOMP         Specify the Vertical Component Channel Identifier.
                            [Default Z].
      --coord-system NAMECONV
                            Coordinate system specification of instrument. (0)
                            Attempt Autodetect between 1 and 2; (1) HZ, HN, HE;
                            (2) Left Handed: HZ, H2 90 CW H1; (3) Right Handed:
                            HZ, H2 90 CCW H1 (4) Left Handed Numeric: H3, H2 90 CW
                            H1. **Note**: This option is not implemented yet. [Default 2]

    Timing Parameters:
      Parameters associated with event timing and window length.

      --start STARTT        Enter Start date for event catalogue search. Note,
                            more recent of this value or station start date will
                            be used.
      --end ENDT            Enter End date for event catalogue search. Note, less
                            recent of this or the station end date will be used.
      --window TWIN         Enter time window length in days. A non-zero value
                            will cause the results to repeat for each set of twin
                            days in the operating window, calculating the change
                            in orientation over time. [Default 0]

    Earthquake Selection Criteria:
      Parameters associated with selecing the subset of earthquakes to use in
      calculations.

      --min-mag MINMAG      Specify the minimum magnitude of Earthquakes to use in
                            the catalogue search. [Default 5.5]
      --min-dist MINDIST    Specify the minimum earthquake distance (in degrees).
                            [Default 5.]
      --max-dist MAXDIST    Specify the maximum earthquake distance (in degrees).
                            [Default 175.]
      --max-dep MAXDEP      Specify maximum Earthquake Depth (km). [Default 150.]

``dl_average``
**************

Description
-----------

Collects the estimated azimuths previously calculated and calculates the
mean value after some quality control thresholding based on the rotated 
waveforms. The error is obtained from a bootstrap analysis of robust estimates.

Usage
-----

.. code-block::

    $ dl_average -h

    #####################################################
    #      _ _                                          #
    #   __| | |    __ ___   _____ _ __ __ _  __ _  ___  #
    #  / _` | |   / _` \ \ / / _ \ '__/ _` |/ _` |/ _ \ #
    # | (_| | |  | (_| |\ V /  __/ | | (_| | (_| |  __/ #
    #  \__,_|_|___\__,_| \_/ \___|_|  \__,_|\__, |\___| #
    #        |_____|                        |___/       #
    #                                                   #
    #####################################################

    usage: dl_average [arguments] <Station Database>

    Program to average the orientations of the seismometer in a station database.

    positional arguments:
      indb                  Station Database to process from.

    optional arguments:
      -h, --help            show this help message and exit
      -V VERB, --verbose VERB
                            Enable Level of verbose output during processing. (0)
                            No Output; (1) Output Event Analysis counter; (2)
                            Counter and results. Default 2
      --load-location LOADLOC
                            Specify Load destination. [Default is DL_RESULTS (and
                            sub-directories based on Station Name)]
      --plot                Plot results at end [Default False]
      --save                Set this option if you wish to save the figure.
                            [Default does not save figure]
      --format FMT          Specify format of figure. Can be any one of the
                            validmatplotlib formats: 'png', 'jpg', 'eps', 'pdf'.
                            [Default 'png']
      --cc CC               Cross-correlation threshold for final estimate.
                            [Default 0.8]
      --min-mag MINMAG      Specify default minimum magnitude to include in average. [Default
                            5.5]
                        
    Station Selection Parameters:
      Parameters to select a specific station.

      --keys STKEYS         Specify list of Station Keys in the database to
                            process.
