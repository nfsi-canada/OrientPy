import numpy as np
from orientpy import arguments
from pkg_resources import resource_filename
from pathlib import Path
from . import get_meta 


dbfile = resource_filename('orientpy',
                           'examples/data/LOBS3.pkl')

def test_dl_calc_args():
    args = arguments.get_dl_calc_arguments([dbfile])
    return args

def test_dl_average_args():
    args = arguments.get_dl_average_arguments([dbfile])
    return args

def test_bng_calc_args():
    args = arguments.get_bng_calc_arguments([dbfile])
    return args

def test_bng_average_args():
    args = arguments.get_bng_average_arguments([dbfile])
    return args

def test_dirs(tmp_path):
    args = test_dl_calc_args()
    db = get_meta.get_stdb()
    stkey = 'YH.LOBS3'

    outdir = tmp_path / args.saveloc 
    if not outdir.exists():
        outdir.mkdir()

    outdir = outdir / stkey.upper()

    if not outdir.exists():
        outdir.mkdir()

