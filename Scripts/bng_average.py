#!/usr/bin/env python

# Final Orientation Calculation File

import stdb
import pickle
import os.path
import numpy as np
from orientpy import utils, arguments, plotting
from pathlib import Path


def main():

    # Run the Input Parser
    args = arguments.get_bng_average_arguments()

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
            filename = indir / folder / "Meta_Data.pkl"
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
                plot.savefig(figname, fmt=args.fmt)
            if args.showplot:
                plot.show()

            plot = plotting.plot_bng_results(stkey, phi, snr, cc, TR, RZ, baz, mag,
                ind, val, err)

            # save figure
            if args.saveplot:
                figname = indir / ('results.' + args.fmt)
                plot.savefig(figname, fmt=args.fmt)
            if args.showplot:
                plot.show()




if __name__ == "__main__":

    # Run main program
    main()
