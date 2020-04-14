#!/usr/bin/env python

# script to test the orient.py running script.

import stdb
import pickle
import os.path
import numpy as np
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from orientpy import BNG, arguments


def main():

    # Run the Input Parser
    args = arguments.get_bng_calc_arguments()

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
        outdir = os.path.join(args.saveloc, stkey.upper())
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

        # Establish client for waveforms
        if len(args.UserAuth) == 0:
            wf_client = Client(args.wf_client)
        else:
            wf_client = Client(args.wf_client,
                               user=args.UserAuth[0],
                               password=args.UserAuth[1])

        # Establish client for catalogue
        cat_client = Client(args.cat_client)

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
                evtdir = outdir + "/" + timekey

                # Check if BNG data already exist and overwrite has been set
                if os.path.isdir(evtdir):
                    if os.path.isfile(evtdir+"/Raw_Data.pkl"):
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
                if not os.path.isdir(evtdir):
                    os.makedirs(evtdir)

                # Save raw Traces
                pickle.dump(bngdata.data, open(
                    evtdir + "/Raw_Data.pkl", "wb"))

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
                pickle.dump(bngdata.meta, open(
                    evtdir + "/Meta_Data.pkl", "wb"))


if __name__ == "__main__":

    # Run main program
    main()
