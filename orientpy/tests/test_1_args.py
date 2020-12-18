import numpy as np
from pkg_resources import resource_filename
from pathlib import Path
import pytest


dbfile = resource_filename('orientpy',
                           'examples/data/LOBS3.pkl')

def test_dl_calc_args():
    from orientpy.scripts import dl_calc as dl
    # no stdb
    with pytest.raises(SystemExit):
        assert dl.get_dl_calc_arguments()

    # defaults
    args0 = dl.get_dl_calc_arguments([dbfile])
    # keys
    args = dl.get_dl_calc_arguments(
        [dbfile, '--keys', 'LOBS3'])
    # start time
    args = dl.get_dl_calc_arguments(
        [dbfile, '--start', '2014-10-01'])
    with pytest.raises(SystemExit):
        assert dl.get_dl_calc_arguments([
            dbfile, '--start', 'abcd'])
    # end time
    args = dl.get_dl_calc_arguments([
        dbfile, '--end', '2015-01-01'])
    with pytest.raises(SystemExit):
        assert dl.get_dl_calc_arguments([
            dbfile, '--end', 'abcd'])
    # user auth.
    args = dl.get_dl_calc_arguments([
        dbfile, '-U', 'bla:bla'])
    with pytest.raises(SystemExit):
        assert dl.get_dl_calc_arguments([
            dbfile, '-U', 'abcd'])
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
    with pytest.raises(SystemExit):
        assert dl.get_dl_average_arguments()
    # defaults
    args0 = dl.get_dl_average_arguments([dbfile])
    # keys
    args = dl.get_dl_average_arguments(
        [dbfile, '--keys', 'LOBS3'])

    return args0

def test_bng_calc_args():
    from orientpy.scripts import bng_calc_auto as bng
    # no stdb
    with pytest.raises(SystemExit):
        assert bng.get_bng_calc_arguments()
    # default
    args0 = bng.get_bng_calc_arguments([dbfile])
    # keys
    args = bng.get_bng_calc_arguments(
        [dbfile, '--keys', 'LOBS3'])
    # start time
    args = bng.get_bng_calc_arguments(
        [dbfile, '--start', '2014-10-01'])
    with pytest.raises(SystemExit):
        assert bng.get_bng_calc_arguments([
            dbfile, '--start', 'abcd'])
    # end time
    args = bng.get_bng_calc_arguments([
        dbfile, '--end', '2015-01-01'])
    with pytest.raises(SystemExit):
        assert bng.get_bng_calc_arguments([
            dbfile, '--end', 'abcd'])
    # user auth.
    args = bng.get_bng_calc_arguments([
        dbfile, '-U', 'bla:bla'])
    with pytest.raises(SystemExit):
        assert bng.get_bng_calc_arguments([
            dbfile, '-U', 'abcd'])
    # local data
    args = bng.get_bng_calc_arguments([
        dbfile, '--local-data', 'bla'])
    # ndval
    args = bng.get_bng_calc_arguments([
        dbfile, '--no-data-zero'])
    # bp
    args = bng.get_bng_calc_arguments([
        dbfile, '--bp', '1.,2.'])
    with pytest.raises(SystemExit):
        assert bng.get_bng_calc_arguments([
            dbfile, '--bp', '1.'])
    # tt
    args = bng.get_bng_calc_arguments([
        dbfile, '--times', '1.,2.'])
    with pytest.raises(SystemExit):
        assert bng.get_bng_calc_arguments([
            dbfile, '--times', '1.'])

    return args0

def test_bng_average_args():
    from orientpy.scripts import bng_average as bng
    # no stdb
    with pytest.raises(SystemExit):
        assert bng.get_bng_average_arguments()
    # defaults
    args0 = bng.get_bng_average_arguments([dbfile])
    # keys
    args = bng.get_bng_average_arguments(
        [dbfile, '--keys', 'LOBS3'])

    return args0
