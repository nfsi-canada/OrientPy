#!/usr/bin/env python

# Final Orientation Calculation File

import stdb
import pickle
import os.path
import numpy as np
from orientpy import utils, arguments, plotting


def main():

    # Run the Input Parser
    args = arguments.get_dl_average_arguments()

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

        # Input directory
        indir = os.path.join(args.loadloc, stkey.upper()) + '/'
        if not os.path.isdir(indir):
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
            filename = indir+"/"+folder+"/Meta_Data.pkl"
            if not os.path.isfile(filename):
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

        if args.showplot:
            plotting.plot_dl_results(stkey, R1phi, R1cc, R2phi, R2cc, ind, 
                val, err, phi, cc, args.cc, loc=indir)


if __name__ == "__main__":

    # Run main program
    main()
