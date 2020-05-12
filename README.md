
![](./orientpy/examples/picture/OrientPy_logo.png)

## Seismic station orientation tools 

OrientPy is a toolbox to help determine seismometer orientation using automated (and manual) 
processing of earthquake data. These methods are particularly useful for broadband 
ocean-bottom seismic stations, but are also applicable to broadband land stations
or shorter period instruments (depending on the method selected). The code uses the 
``StDb`` package for querying and building a station database and can be used through 
command-line scripts.Currently the toolbox includes the following methods: 

- **DL** (Doran and Laske, 2017): Based on Rayleigh-wave polarization at a range of
  frequencies and for the two fundamental mode Rayleigh wave orbits. 

- **BNG** (Braunmiller, Nabelek and Ghods, 2020): Based on P-wave polarization from 
  regional and teleseismic earthquakes.  

- **LKSS** (Lim et al., 2018): Based on the harmonic decomposition of radial and 
  transverse receiver functions near zero lag times.

Each method can be used independently to produce an estimate of station orientation, in
terms of the azimuth of seismic component `1` (or `N`).

[![Build Status](https://travis-ci.com/nfsi-canada/OrientPy.svg?branch=master)](https://travis-ci.com/nfsi-canada/OrientPy)
[![codecov](https://codecov.io/gh/nfsi-canada/OrientPy/branch/master/graph/badge.svg)](https://codecov.io/gh/nfsi-canada/OrientPy)

Installation, Usage, API documentation and scripts are described at 
https://nfsi-canada.github.io/OrientPy/.

#### Contributing

All constructive contributions are welcome, e.g. bug reports, discussions or suggestions for new features. You can either [open an issue on GitHub](https://github.com/nfsi-canada/OrientPy/issues) or make a pull request with your proposed changes. Before making a pull request, check if there is a corresponding issue opened and reference it in the pull request. If there isn't one, it is recommended to open one with your rationale for the change. New functionality or significant changes to the code that alter its behavior should come with corresponding tests and documentation. If you are new to contributing, you can open a work-in-progress pull request and have it iteratively reviewed.

Examples of straightforward contributions include editing the documentation or adding notebooks/scripts that describe example usage of the code in publications. Suggestions for improvements (speed, accuracy, flexibility, etc.) are also welcome.

#### References

- Braunmiller, J., Nabelek, J., and Ghods, A., 2020, Sensor orientation of Iranian broadband
  seismic stations from P-wave particle motion, *Seism. Res. Lett.*, doi:10.1785/0220200019.

- Doran, A. K., and Laske, G., 2017, Ocean-bottom seismometer instrument orientation 
  via automated Rayleigh-wave arrival-angle measurements, *Bull. Seism. Soc. Am.*,
  107, 691-708.

- Lim, H., Kim, Y., Song, T.-R. A., and Shen, Z., 2018, Measurement of 
  seismometer orientation using the tangential P-wave
  receiver function based on harmonic decomposition, *Geophys. J. Int.*, 212,
  1747-1765.

