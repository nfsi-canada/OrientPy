
.. figure:: ../orientpy/examples/picture/OrientPy_logo.png
   :align: center

Licence
=======

Copyright 2020 Pascal Audet 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Installation
============

Dependencies
------------

The current version has been tested using **Python3.6 and 3.7** \
Also, the following packages are required:

- `stdb <https://github.com/schaefferaj/StDb>`_
- `geographiclib <https://geographiclib.sourceforge.io/html/python/>`_

Other required packages (e.g., ``obspy``)
will be automatically installed by ``stdb``.

Conda environment
-----------------

We recommend creating a custom ``conda`` environment
where ``OrientPy`` can be installed along with some of its dependencies.

.. sourcecode:: bash

   conda create -n orient python=3.7 obspy -c conda-forge

Activate the newly created environment:

.. sourcecode:: bash

   conda activate orient

Install remaining dependencies using ``pip`` inside the ``orient`` environment:

.. sourcecode:: bash

   pip install stdb
   pip install geographiclib


Installing from Pypi
--------------------

*This option is not available at this time*

Installing from source
----------------------

- Clone the repository:

.. sourcecode:: bash

   git clone https://github.com/paudetseis/OrientPy.git
   cd OrientPy

- Install using pip:

.. sourcecode:: bash

   pip install .

.. Basic Usage
.. ===========

.. Calculating Receiver Functions
.. ------------------------------

.. The basic class of ``RfPy`` is :class:`~rfpy.rfdata.RFData`, which contains attributes and 
.. methods for the calculation of single-station, teleseismic 
.. `P`-wave receiver functions from three-component seismograms. A :class:`~rfpy.rfdata.RFData`
.. object contains three main attributes: a :class:`~stdb.StDb` object with station information,
.. a :class:`~rfpy.rfdata.Meta` object containing event meta data, and a :class:`~obspy.core.Stream`
.. object containing the unrotated 3-component seismograms. Additional processing attributes 
.. are added as the analysis proceeds. The sequence of initialization and addition of attributes 
.. is important, as described in the documentation below. 

.. Note that, at the end of the process, the :class:`~rfpy.rfdata.RFData` object will further contain
.. a :class:`~obspy.core.Stream` object as an additional attribute, containing the receiver function data.

.. .. note::

..     A :class:`~rfpy.rfdata.RFData` object is meant to facilitate processing of single-station 
..     and single-event P-wave receiver functions. For processing multiple event-station pairs, 
..     an equal number of :class:`~rfpy.rfdata.RFData` objects need to be 
..     created. See the accompanying Scripts for details.

.. Initialization
.. ++++++++++++++

.. A ``RFData`` object is initialized with an :class:`~stdb.StDb` object, e.g. consider such an 
.. object ``sta``:

.. .. sourcecode:: python

..     >>> from rfpy import RFData
..     >>> rfdata = RFData(sta)


.. Once the object is initialized, the first step is to add an :class:`obspy.core.event.Event``
.. object. For example, given such an object ``ev``:

.. .. sourcecode:: python

..     >>> rfdata.add_event(ev)

.. Now that the event has been added, the :class:`~rfpy.rfdata.RFData` object has determined
.. whether or not it is suitable for receiver function analysis (i.e., 
.. if the event is within a suitable epicentral distance range), which is
.. available as a new ``meta`` attribute:

.. .. sourcecode:: python

..     >>> rfdata.meta.accept
..     True

.. .. note::

..     Alternatively, the :func:`~rfpy.rfdata.RFData.add_event` 
..     (or :func:`~rfpy.rfdata.RFData.add_data`) method can be used
..     with the argument ``returned=True`` to return the ``accept`` attribute
..     directly.

..     .. sourcecode:: python

..         >>> rfdata.add_event(ev, returned=True)
..         True

.. If the ``accept`` attribute is ``True``, continue with the analysis by
.. adding raw three-component data. There are two methods to perform this step.
.. If the data are available in memory (e.g., in a :class:`~obspy.core.Stream` object ``stream``), 
.. one can use the :func:`~rfpy.rfdata.RFData.add_data` method directly:

.. .. sourcecode:: python

..     >>> rfdata.add_data(stream)

.. .. warning::

..     **Do not** simply add a :class:`~obspy.core.Stream` object as an 
..     attribute ``data`` to the :class:`~rfpy.rfdata.RFData`
..     object (e.g., ``rfdata.data = stream``). Instead use this method, as it checks 
..     whether or not the data are suitable for receiver function analysis.

.. Otherwise, one can use the method :func:`~rfpy.rfdata.RFData.download_data` to obtain 
.. the three-component data from an FDSN Client: 

.. .. sourcecode:: python

..     >>> rfdata.download_data(client)

.. The ``accept`` attribute will be updated with the availability of the ``data``
.. attribute, i.e. if no data is available, the ``accept`` attribute will be set
.. to ``False``. The methods to add data can also be used with the argument 
.. ``returned=True`` to report whether or not the data are available. 

.. Receiver function processing
.. ++++++++++++++++++++++++++++

.. Now that we have complete meta data and raw seismogram data, we can
.. use methods to rotate and/or calculate the signal-to-noise ratio. 
.. The rotation flag is set in the ``rfdata.meta.align`` attribute, which by
.. default is ``'ZRT'``. This means that ``'ZNE'`` data will be rotated to ``'ZRT'`` 
.. before deconvolution, automatically. However, we can set a different alignment
.. (e.g., ``'LQT'`` or ``'PVH'``) and perform the rotation prior to deconvolution.
.. Once rotation is performed, however, the initial ``'ZNE'`` data is no longer 
.. available and further rotation cannot be performed:

.. .. sourcecode:: python

..     >>> rfdata.rotate()        
..     >>> rfdata.meta.rotated
..     True
..     >>> rfdata.meta.align
..     'ZRT'
..     >>> rfdata.rotate(align='PVH')
..     ...
..     Exception: Data have been rotated already - aborting

.. The SNR is calculated based on the ``align`` attribute, on the first component
.. (e.g., either ``'Z'``, ``'L'`` or ``'P'``). Therefore, this method is typically
.. carried out following the ``rotate`` method:

.. .. sourcecode:: python

..     >>> rfdata.calc_snr()
..     >>> type(rfdata.meta.snr)
..     float

.. Finally, the last step is to perform the deconvolution using the method 
.. :func:`~rfpy.rfdata.RFData.deconvolve`,
.. which stores the receiver function data as a new attribute ``rf``, which is a 
.. three-component :class:`~obspy.core.Stream` object:

.. .. sourcecode:: python

..     >>> rfdata.deconvolve()

.. Although no plotting method is provided for the :class:`~rfpy.rfdata.RFData` object, 
.. the ``rf`` attribute is a :class:`~obspy.core.Stream`
.. object that can be plotted using the :func:`~rfpy.rfdata.RFData.plot` method 
.. (e.g., ``rfdata.rf.plot()``).

.. Following receiver function deconvolution, all the information is stored in the attributes 
.. of the object. Ultimately, a method is available to convert the ``RFData`` object to a
.. :class:`~obspy.core.Stream` object with new attributes:

.. .. sourcecode:: python

..     >>> rfstream = rfdata.to_stream()

.. Demo example
.. ++++++++++++

.. To look at a concrete example for station `MMPY <http://ds.iris.edu/mda/NY/MMPY/>`_, 
.. consider the demo data provided with the package and process them using all default values: 

.. .. sourcecode:: python

..     >>> from rfpy import RFData
..     >>> rfdata = RFData('demo')
..     Uploading demo station data - station NY.MMPY

.. Check out its attributes (initialization only stores the ``sta`` attribute)

.. .. sourcecode:: python

..     >>> rfdata.__dict__
..     {'sta': {'station': 'MMPY',
..       'network': 'NY',
..       'altnet': [],
..       'channel': 'HH',
..       'location': ['--'],
..       'latitude': 62.618919,
..       'longitude': -131.262466,
..       'elevation': 0.0,
..       'startdate': 2013-07-01T00:00:00.000000Z,
..       'enddate': 2599-12-31T23:59:59.000000Z,
..       'polarity': 1.0,
..       'azcorr': 0.0,
..       'status': 'open'},
..      'meta': None,
..      'data': None}

.. Now import an event:

.. .. sourcecode:: python

..     >>> rfdata.add_event('demo')
..     2015-07-03T06:43:22.840000Z | +10.139, +125.977 | 6.1 MW

.. Print the content of the object meta data

.. .. sourcecode:: python

..     >>> rfdata.meta.__dict__
..   {'time': 2015-07-03T06:43:22.840000Z,
..    'lon': 125.977,
..    'lat': 10.1385,
..    'dep': 50.8,
..    'mag': 6.1,
..    'epi_dist': 9654.95322183263,
..    'az': 26.754776375134075,
..    'baz': 286.03267562515674,
..    'gac': 86.82908036528741,
..    'ttime': 766.29205941718612,
..    'ph': 'P',
..    'slow': 0.043863956666625202,
..    'inc': 14.712262941689376,
..    'accept': True,
..    'vp': 6.0,
..    'vs': 3.6,
..    'align': 'ZRT',
..    'rotated': False,
..    'snr': None}

.. .. note::

..     Once the event object is loaded, it is possible to edit the attributes
..     of ``meta``, although we recommend only editing ``vp``, ``vs`` or 
..     ``align``, and avoid editing any of the station-event attributes

..     .. sourcecode:: python

..         >>> rfdata.meta.vp = 5.5
..         >>> rfdata.meta.vs = 3.3
..         >>> rfdata.meta.vp, rfdata.meta.vs
..         (5.5, 3.3)
..         >>> rfdata.meta.align = 'LQT'
..         >>> rfdata.meta.align
..         'LQT'

.. Now add data to the object:

.. .. sourcecode:: python

..     >>> rfdata.add_data('demo')
..   3 Trace(s) in Stream:
..   NY.MMPY..HHZ | 2015-07-03T06:54:09.140000Z - 2015-07-03T06:58:08.940000Z | 5.0 Hz, 1200 samples
..   NY.MMPY..HHN | 2015-07-03T06:54:09.140000Z - 2015-07-03T06:58:08.940000Z | 5.0 Hz, 1200 samples
..   NY.MMPY..HHE | 2015-07-03T06:54:09.140000Z - 2015-07-03T06:58:08.940000Z | 5.0 Hz, 1200 samples

.. Perform receiver function deconvolution using default values:

.. .. sourcecode:: python

..     >>> rfdata.deconvolve()
..     Warning: Data have not been rotated yet - rotating now
..     Warning: SNR has not been calculated - calculating now using 

..     >>> rfdata.rf
..     3 Trace(s) in Stream:
..     NY.MMPY..RFZ | 2015-07-03T06:56:04.140000Z - 2015-07-03T06:57:59.140000Z | 5.0 Hz, 576 samples
..     NY.MMPY..RFR | 2015-07-03T06:56:04.140000Z - 2015-07-03T06:57:59.140000Z | 5.0 Hz, 576 samples
..     NY.MMPY..RFT | 2015-07-03T06:56:04.140000Z - 2015-07-03T06:57:59.140000Z | 5.0 Hz, 576 samples    

..     >>> rfstream = rfdata.to_stream()
..     >>> rfstream
..     3 Trace(s) in Stream:
..     NY.MMPY..RFZ | 2015-07-03T06:56:04.140000Z - 2015-07-03T06:57:59.140000Z | 5.0 Hz, 576 samples
..     NY.MMPY..RFR | 2015-07-03T06:56:04.140000Z - 2015-07-03T06:57:59.140000Z | 5.0 Hz, 576 samples
..     NY.MMPY..RFT | 2015-07-03T06:56:04.140000Z - 2015-07-03T06:57:59.140000Z | 5.0 Hz, 576 samples

.. Check out new stats in traces

.. .. sourcecode:: python

..     >>> rfstream[0].stats.snr
..     10.799614447256117
..     >>> rfstream[0].stats.slow
..     0.043863956666625202
..     >>> rfstream[0].stats.baz
..     286.03267562515674
..     >>> rfstream[0].stats.is_rf
..     True

.. Plot filtered and trimmed ``rfstream``


.. ..sourcecode:: python

..     >>> rfstream.filter('bandpass', freqmin=0.05, freqmax=0.5)
..     >>> t1 = rfstream[0].stats.starttime
..     >>> rfstream.trim(t1, t1+30.)
..     >>> rfstream.plot()

.. .. figure:: ../rfpy/examples/data/Figure_rfdata_demo.png
..    :align: center


.. Post-Processing: `H-k` stacking
.. -------------------------------

.. The class :class:`~rfpy.hk.HkStack` contains attributes and methods to calculate thickness (`H`) 
.. and Vp/Vs ratio (`k`) of the crust (in reality, `H` refers to Moho depth, and `k` is Vp/Vs of 
.. the medium from the surface to `H`) based on moveout times of direct `Ps` and reverberated 
.. `Pps` and `Pss` phases from radial-component receiver functions. The individual 
.. phase stacks are obtained from the median weighted by the phase of individual 
.. signals. Methods are available to combine the phase stacks into a weighted sum
.. or a product. 

.. Initialization
.. ++++++++++++++

.. A ``HkStack`` object is initialized with a :class:`~obspy.core.Stream` 
.. object containing radial receiver function data. The :class:`~obspy.core.Stream` 
.. is built by adding (or appending) radial receiver functions obtained from valid
.. :class:`~rfpy.rfdata.RFData` objects using the :func:`~rfpy.rfdata.RFData.to_stream`
.. method.

.. .. sourcecode:: python

..     >>> from rfpy import HkStack
..     >>> hkstack = HkStack(rfstream)

.. The ``rfstream`` typically requires minimal pre-processing, such as
.. bandpass filtering to enhance the converted and reverberated phases.
.. For example:

.. .. sourcecode:: python

..     >>> rfstream.filter('bandpass', freqmin=0.05, freqmax=0.75, corners=2, zerophase=True)
..     >>> hkstack = HkStack(rfstream)

.. .. note::

..     It is also possible to use two ``rfstream`` objects during initialization
..     of the :class:`~rfpy.hk.HkStack` object - one for the direct conversion 
..     (i.e., ``'ps'`` phase), 
..     and the second one for the reverberated phases (i.e., ``'pps'``, ``'pss'``).
..     The second ``rfstream`` should therefore be a copy of the first one, but perhaps
..     filtered uding different frequency corners:

..     .. sourcecode:: python

..         >>> rfstream2 = rfstream.copy()
..         >>> rfstream2.filter('bandpass', freqmin=0.05, freqmax=0.35, corners=2, zerophase=True)
..         >>> hkstack = HkStack(rfstream, rfstream2)

.. To speed things up during processing (and to avoid redundant stacking), it is possible to
.. use one of the :func:`~rfpy.binning` functions, alghouth **not** the 
.. :func:`~rfpy.binning.bin_all` function, e.g.,

.. .. sourcecode:: python

..     >>> from rfpy.binning import bin
..     >>> rfstream_binned = rfstream.bin(typ='slow', nbin=21)
..     >>> hkstack = HkStack(rfstream_binned)

.. H-k processing
.. ++++++++++++++

.. Once the :class:`~rfpy.hk.HkStack` object is initialized with the ``rfstream``, a findividual phase
.. stacks can be calculated automatically using the default settings:

.. .. sourcecode:: python

..     >>> hkstack.stack()

.. The only parameter to set is the `P`-wave velocity of the crust - if not set,
.. the default value of 6.0 km/s is used (available as the attribute ``hkstack.vp``).
.. To change the search bounds for the phase stacks, we can edit the attributes of the
.. :class:`~rfpy.hk.HkStack` object prior to calling the method :func:`~rfpy.hk.HkStack.stack`:

.. .. sourcecode:: python

..     >>> hkstack.hbound = [15., 40.]
..     >>> hkstack.dh = 1.5
..     >>> hkstack.kbound = [1.6, 2.0]
..     >>> hkstack.dk = 0.01
..     >>> hkstack.stack(vp=5.5)

.. .. warning::

..     Setting small values for ``hkstack.dh`` and ``hkstack.dk`` will slow down
..     the processing significantly, but produce much cleaner and more precise
..     stacks.

.. In the presence of a dipping Moho interface, it is possible to use the method
.. :func:`~rfpy.hk.HkStack.stack_dip`, with the additional ``strike`` and ``dip`` arguments.
.. If not specified, the code will use the default values stored as attributes of the
.. :class:`~rfpy.hk.HkStack` object:

.. .. sourcecode:: python

..     >>> hkstack.stack_dip(strike=215., dip=25., vp=5.5)

.. Once the phase stacks are calculated and stored as attributes of the object,
.. we can call the method :func:`~rfpy.hk.HkStack.average` to combine the phase stacks
.. into a single, final stack. By default the final stack is a simple weighted sum 
.. of the individual phase stacks, using weights defined as object attributes:

.. .. sourcecode:: python

..     >>> hkstack.weights
..     [0.5, 2., -1.]
..     >>> hkstack.average()

.. To produce a final stack that consists of the product of the positive parts
.. of individual phase stacks (to enhance normal-polarity Moho arrivals and ignore
.. un-modelled negative polarity signals), use the ``typ='prod'`` argument:

.. .. sourcecode:: python

..     >>> hkstack.average(typ='prod')

.. The estimates of `H` and `k` are determined from the maximum value in the final
.. stack as attributes ``hkstack.h0`` and ``hkstack.k0``. The method will also 
.. call the :func:`~rfpy.hk.HkStack.error` method to calculate the errors
.. and error contour around the solution.

.. The individual and final stacks can be plotted by calling the method 
.. :func:`~rfpy.hk.HkStack.plot`:

.. .. sourcecode:: python

..     >>> hkstack.plot()


.. Demo example
.. ++++++++++++

.. Initialize object with demo data for station `MMPY <http://ds.iris.edu/mda/NY/MMPY/>`_:

.. .. sourcecode:: python

..     >>> from rfpy import HkStack
..     >>> hkstack = HkStack('demo')
..     Uploading demo data - station NY.MMPY

..     >>> # Check content of object
..     >>> hkstack.__dict__
..     {'rfV1': 66 Trace(s) in Stream:

..     NY.MMPY..RFV | 2016-05-31T10:11:49.520000Z - 2016-05-31T10:13:44.520000Z | 5.0 Hz, 576 samples
..     ...
..     (64 other traces)
..     ...
..     NY.MMPY..RFV | 2015-06-08T06:10:13.330000Z - 2015-06-08T06:12:08.330000Z | 5.0 Hz, 576 samples

..     [Use "print(Stream.__str__(extended=True))" to print all Traces],
..      'rfV2': None,
..      'strike': None,
..      'dip': None,
..      'vp': 6.0,
..      'kbound': [1.56, 2.1],
..      'dk': 0.02,
..      'hbound': [20.0, 50.0],
..      'dh': 0.5,
..      'weights': [0.5, 2.0, -1.0],
..      'phases': ['ps', 'pps', 'pss']}

.. These receiver functions have been obtained by adding :class:`~rfpy.rfdata.RFData` objects
.. as streams to an :class:`~obspy.core.Stream` object, without other processing. Note that they
.. are aligned in the ``'PVH'`` coordinate system, as specified in the channel name (i.e., ``'RFV'`` for
.. the radial component). To prepare them for stacking, we can bin the receiver functions into
.. back-azimuth and slowness bins (in the presence of a dipping interface), or simply slowness bins 
.. (for horizontal interfaces):

.. .. sourcecode:: python

..     >>> from rfpy import binning
..     >>> rfV_binned = binning.bin(hkstack.rfV1, typ='slow', nbin=21)[0]
..     >>> hkstack.rfV1 = rfV_binned

.. it is straightforward to directly
.. filter the :class:`~obspy.core.Stream` object, and perhaps also add a copy of the stream
.. with a different frequency corner as another attribute ``rfV2``, as suggested above:

.. .. sourcecode:: python

..     >>> hkstack.rfV2 =  hkstack.rfV1.copy()
..     >>> hkstack.rfV1.filter('bandpass', freqmin=0.05, freqmax=0.5, corners=2, zerophase=True)
..     >>> hkstack.rfV2.filter('bandpass', freqmin=0.05, freqmax=0.35, corners=2, zerophase=True)

.. Now simply process the hkstack object using the default values to obtain `H` and `k` estimates

.. .. sourcecode:: python

..     >>> hkstack.stack()
..     Computing: [###############] 61/61

..     >>> hkstack.average()

.. The final estimates are available as attributes

.. .. sourcecode:: python

..     >>> hkstack.h0
..     34.0
..     >>> hkstack.err_h0
..     5.5
..     >>> hkstack.k0
..     1.72
..     >>> hkstack.err_k0
..     0.14

.. Plot the stacks with error contours

.. .. sourcecode:: python

..     >>> hkstack.plot()

.. .. figure:: ../rfpy/examples/data/Figure_hk_demo.png
..    :align: center


.. Post-Processing: Harmonic Decomposition
.. ---------------------------------------

.. The class :class:`~rfpy.harmonics.Harmonics` contains attributes and methods to 
.. calculate the first five 
.. harmonic components of radial and transverse component receiver function
.. data from a singular value decomposition. The harmonic decomposition can 
.. be performed at a fixed azimuth (i.e., along some known dominant strike 
.. direction in the subsurface), or alternatively the decomposition can 
.. be optimized to search for the dominant azimuth that maximizes the energy
.. on one of the components. This direction can be interpreted as the 
.. strike of a dipping interface or can be related to anisotropic axes.

.. Initialization
.. ++++++++++++++

.. a :class:`~rfpy.harmonics.Harmonics` object is initialized with **both** radial
.. and transverse component receiver function :class:`~obspy.core.Stream` objects.
.. The :class:`~obspy.core.Stream` objects are built by adding (or appending) 
.. radial and transverse receiver functions obtained from valid
.. :class:`~rfpy.rfdata.RFData` objects using the :func:`~rfpy.rfdata.RFData.to_stream`
.. method.

.. .. sourcecode:: python

..     >>> from rfpy import Harmonics
..     >>> harmonics = Harmonics(rfRstream, rfTstream)

.. .. note::

..     The ``rfRstream`` and ``rfTstream`` typically require minimal pre-processing, such as
..     bandpass filtering to enhance the converted and reverberated phases.
..     For example:

..     .. sourcecode:: python

..         >>> rfRstream.filter('bandpass', freqmin=0.05, freqmax=0.75, corners=2, zerophase=True)
..         >>> rfTstream.filter('bandpass', freqmin=0.05, freqmax=0.75, corners=2, zerophase=True)
..         >>> harmonics = Harmonics(rfRstream, rfTstream)

.. .. warning::

..     The radial and transverse components should not be mixed, and should contain 
..     purely radial and purely transverse components (i.e. no mixing of components). 
..     Furthermore, the :class:`~obspy.core.Stream` objects should have equal length
..     and the same ordering.

.. Harmonic decomposition
.. ++++++++++++++++++++++

.. Once the :class:`~rfpy.harmonics.Harmonics` object is initialized, processing is done by typing:

.. .. sourcecode:: python

..     >>> harmonics.dcomp_fix_azim() 

.. Or, alternatively,

.. .. sourcecode:: python

..     >>> harmonics.dcomp_find_azim()

.. In either case the harmonic components are available as an attribute of type
.. :class:`~obspy.core.Stream` (``harmonics.hstream``) and, if available, the azimuth
.. of the dominant direction (``harmonics.azim``). 

.. .. note::

..     When using the method :func:`rfpy.harmonics.dcomp_find_azim`, it is possible to
..     specify a range of values over which to perform the search using the arguments
..     ``xmin`` and ``xmax``, where `x` refers to the independent variable (i.e., time
..     or depth, if the streams have been converted from time to depth a priori). 

.. Once the harmonic decomposition is performed, the components can be plotted using
.. the method :func:`~rfpy.harmonics.Harmonics.plot`

.. .. sourcecode:: python

..     >>> harmonics.plot()

.. Forward modeling
.. ++++++++++++++++

.. If the ``hstream`` attribute is available, it is possible to *forward model* receiver functions
.. for a range of back-azimuth values, or just a single value. In case the back-azimuths are
.. not specified, the method will use the range of values available in the original
.. radial and transverse component receiver function data.

.. .. sourcecode:: python

..     >>> harmonics.forward()

.. The new `predicted` radial and transverse component receiver functions are available
.. as attributes of type :class:`~obspy.core.Stream` (``harmonics.forwardR`` and ``harmonics.forwardT``)

.. Demo example
.. ++++++++++++

.. Initialize object with demo data for station `MMPY <http://ds.iris.edu/mda/NY/MMPY/>`_:

.. .. sourcecode:: python

..     >>> from rfpy import Harmonics
..     >>> harmonics = Harmonics('demo')
..     Uploading demo data - station NY.MMPY

..     >>> # Check content of object
..     >>> harmonics.__dict__
..     {'strV': 66 Trace(s) in Stream:

..     NY.MMPY..RFV | 2016-05-31T10:11:49.520000Z - 2016-05-31T10:13:44.520000Z | 5.0 Hz, 576 samples
..     ...
..     (64 other traces)
..     ...
..     NY.MMPY..RFV | 2015-06-08T06:10:13.330000Z - 2015-06-08T06:12:08.330000Z | 5.0 Hz, 576 samples

..     [Use "print(Stream.__str__(extended=True))" to print all Traces],
..      'strH': 66 Trace(s) in Stream:

..     NY.MMPY..RFH | 2016-05-31T10:11:49.520000Z - 2016-05-31T10:13:44.520000Z | 5.0 Hz, 576 samples
..     ...
..     (64 other traces)
..     ...
..     NY.MMPY..RFH | 2015-06-08T06:10:13.330000Z - 2015-06-08T06:12:08.330000Z | 5.0 Hz, 576 samples

..     [Use "print(Stream.__str__(extended=True))" to print all Traces],
..      'azim': 0,
..      'xmin': 0.0,
..      'xmax': 40.0}

.. As with the :class:`~rfpy.hk.HkStack` object, these receiver functions have been obtained 
.. by adding :class:`~rfpy.rfdata.RFData` objects
.. as streams to an :class:`~obspy.core.Stream` object, without other processing. Note that they
.. are aligned in the ``'PVH'`` coordinate system, as specified in the channel name (i.e., ``'RFV'`` 
.. and ``'RFH'``). To prepare them for harmonic decomposition, we can bin the receiver functions into
.. back-azimuth and slowness bins :

.. .. sourcecode:: python

..     >>> from rfpy import binning
..     >>> str_binned = binning.bin_baz_slow(harmonics.radialRF, harmonics.transvRF)
..     >>> harmonics.radialRF = str_binned[0]
..     >>> harmonics.transvRF = str_binned[1]

.. It is straightforward to directly
.. filter the :class:`~obspy.core.Stream` object, and perhaps also add a copy of the stream
.. with a different frequency corner as another attribute ``rfV2``, as suggested above:

.. .. sourcecode:: python

..     >>> harmonics.radialRF.filter('bandpass', freqmin=0.05, freqmax=0.5, corners=2, zerophase=True)
..     >>> harmonics.transvRF.filter('bandpass', freqmin=0.05, freqmax=0.5, corners=2, zerophase=True)

.. Now simply perform harmonic decomposition

.. .. sourcecode:: python

..     >>> harmonics.dcomp_fix_azim()
..     Decomposing receiver functions into baz harmonics for azimuth =  0

.. Plot them

.. .. sourcecode:: python

..     >>> harmonics.plot(ymax=10.)

.. .. figure:: ../rfpy/examples/data/Figure_harmonics_demo.png
..    :align: center
