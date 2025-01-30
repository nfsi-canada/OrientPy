import numpy as np
from pkg_resources import resource_filename
from pathlib import Path
import pytest


dbfile = resource_filename('orientpy',
                           'examples/data/LOBS3.pkl')

def test_dl_calc_args():
    from orientpy.scripts import dl_calc as dl
    # no stdb
    with pytest.raises(SystemExit) as excinfo:
        dl.get_dl_calc_arguments()
        assert excinfo.value.code == 1

    # defaults
    args0 = dl.get_dl_calc_arguments([dbfile])
    # keys
    args = dl.get_dl_calc_arguments(
        [dbfile, '--keys', 'LOBS3'])
    # start time
    args = dl.get_dl_calc_arguments(
        [dbfile, '--start', '2014-10-01'])
    with pytest.raises(SystemExit) as excinfo:
        dl.get_dl_calc_arguments([
            dbfile, '--start', 'abcd'])
        assert excinfo.value.code == 1
    # end time
    args = dl.get_dl_calc_arguments([
        dbfile, '--end', '2015-01-01'])
    with pytest.raises(SystemExit) as excinfo:
        dl.get_dl_calc_arguments([
            dbfile, '--end', 'abcd'])
        assert excinfo.value.code == 1
    # user auth.
    args = dl.get_dl_calc_arguments([
        dbfile, '-U', 'bla:bla'])
    with pytest.raises(SystemExit) as excinfo:
        dl.get_dl_calc_arguments([
            dbfile, '-U', 'abcd'])
        assert excinfo.value.code == 1
    # local data
    args = dl.get_dl_calc_arguments([
        dbfile, '--local-data', 'bla'])
    # ndval
    args = dl.get_dl_calc_arguments([
        dbfile, '--no-data-zero'])

    return args0

def test_dl_average_args():
    from orientpy.scripts import dl_average as dl
    # no stdb
    with pytest.raises(SystemExit) as excinfo:
        dl.get_dl_average_arguments()
        assert excinfo.value.code == 1
    # defaults
    args0 = dl.get_dl_average_arguments([dbfile])
    # keys
    args = dl.get_dl_average_arguments(
        [dbfile, '--keys', 'LOBS3'])

    return args0

def test_bng_calc_args():
    from orientpy.scripts import bng_calc_auto as bng
    # no stdb
    with pytest.raises(SystemExit) as excinfo:
        bng.get_bng_calc_arguments()
        assert excinfo.value.code == 1
    # default
    args0 = bng.get_bng_calc_arguments([dbfile])
    # keys
    args = bng.get_bng_calc_arguments(
        [dbfile, '--keys', 'LOBS3'])
    # start time
    args = bng.get_bng_calc_arguments(
        [dbfile, '--start', '2014-10-01'])
    with pytest.raises(SystemExit) as excinfo:
        bng.get_bng_calc_arguments([
            dbfile, '--start', 'abcd'])
        assert excinfo.value.code == 1
    # end time
    args = bng.get_bng_calc_arguments([
        dbfile, '--end', '2015-01-01'])
    with pytest.raises(SystemExit) as excinfo:
        bng.get_bng_calc_arguments([
            dbfile, '--end', 'abcd'])
        assert excinfo.value.code == 1
    # user auth.
    args = bng.get_bng_calc_arguments([
        dbfile, '-U', 'bla:bla'])
    with pytest.raises(SystemExit) as excinfo:
        bng.get_bng_calc_arguments([
            dbfile, '-U', 'abcd'])
        assert excinfo.value.code == 1
    # local data
    args = bng.get_bng_calc_arguments([
        dbfile, '--local-data', 'bla'])
    # ndval
    args = bng.get_bng_calc_arguments([
        dbfile, '--no-data-zero'])
    # bp
    args = bng.get_bng_calc_arguments([
        dbfile, '--bp', '1.,2.'])
    with pytest.raises(SystemExit) as excinfo:
        bng.get_bng_calc_arguments([
            dbfile, '--bp', '1.'])
        assert excinfo.value.code == 1
    # tt
    args = bng.get_bng_calc_arguments([
        dbfile, '--times', '1.,2.'])
    with pytest.raises(SystemExit) as excinfo:
        bng.get_bng_calc_arguments([
            dbfile, '--times', '1.'])
        assert excinfo.value.code == 1

    return args0

def test_bng_average_args():
    from orientpy.scripts import bng_average as bng
    # no stdb
    with pytest.raises(SystemExit) as excinfo:
        bng.get_bng_average_arguments()
        assert excinfo.value.code == 1
    # defaults
    args0 = bng.get_bng_average_arguments([dbfile])
    # keys
    args = bng.get_bng_average_arguments(
        [dbfile, '--keys', 'LOBS3'])

    return args0
