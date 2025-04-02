
.. figure:: ../orientpy/examples/picture/OrientPy_logo.png
   :align: center

Documentation
=============

``OrientPy`` is a package containing Python tools to estimate seismometer 
orientation using automated (and manual) processing of earthquake data. These 
methods are particularly useful for broadband ocean-bottom seismic stations, 
but are also applicable to broadband land stations or shorter period instruments 
(depending on the method selected). Currently the toolbox includes the following methods: 

- **DL** (Doran and Laske, 2017): Based on Rayleigh-wave polarization at a range of
  frequencies and for the two fundamental mode Rayleigh wave orbits. 

- **BNG** (Braunmiller, Nabelek and Ghods, 2020): Based on P-wave polarization of 
  earthquakes at regional teleseismic distances.  

Each method can be used independently to produce an estimate of station orientation, in
terms of the azimuth of seismic component '1' (or 'N'). The code uses 
the ``StDb`` package for querying and building a station database 
used in command-line scripts. Tutorials are provided in the sections below.

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3905404.svg
   :target: https://doi.org/10.5281/zenodo.3905404
.. image:: https://github.com/nfsi-canada/OrientPy/workflows/Build/badge.svg
   :target: https://github.com/nfsi-canada/OrientPy/actions
.. image:: https://codecov.io/gh/nfsi-canada/OrientPy/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/nfsi-canada/OrientPy

.. note::

   The method ``DL`` is heavily based on the code ``DLOPy`` by
   Doran and Laske (2017). The original code can be found 
   `here <https://igppweb.ucsd.edu/~adoran/DLOPy.html>`_

.. toctree::
   :maxdepth: 1
   :caption: Quick Links

   links

.. toctree::
   :maxdepth: 1
   :caption: Disclaimer

   disclaimer

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   orientpy

.. toctree::
   :maxdepth: 2
   :caption: API Documentation

   api

.. toctree::
   :maxdepth: 2
   :caption: Scripts & Tutorials

   scripts
   tutorials