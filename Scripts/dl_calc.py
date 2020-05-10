#!/usr/bin/env python

# script to test the orient.py running script.

import stdb
import pickle
import os.path
import numpy as np
from obspy.clients.fdsn import Client
from orientpy import DL, arguments, utils
from pathlib import Path


def main():

    # Run the Input Parser
    args = arguments.get_dl_calc_arguments()

    # Load Database
    db = stdb.io.load_db(fname=args.indb)

    # Construct station key loop
    allkeys = db.keys()
    sorted(allkeys)

    # Extract key subset
    if len(args.stkeys) > 0:
        stkeys = []
        for skey in args.stkeys:
            stkeys.extend([s for s in allkeys if skey in s])
    else:
        stkeys = db.keys()
        sorted(stkeys)

    # Loop over station keys
    for stkey in list(stkeys):

        sta = db[stkey]

        # Output directory
        outdir = Path(args.saveloc) / Path(stkey.upper())
        if not outdir.exists():
            outdir.mkdir()

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

                # Create Folder if it doesn't exist
                if not evtdir.exists():
                    evtdir.mkdir()

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
