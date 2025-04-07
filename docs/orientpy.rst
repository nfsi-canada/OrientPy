
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

The current version has been tested using **Python3.10** \
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

   conda create -n orient "python=3.10" "setuptools=60" obspy -c conda-forge

Activate the newly created environment:

.. sourcecode:: bash

   conda activate orient

Install remaining dependencies using ``pip`` inside the ``orient`` environment:

.. sourcecode:: bash

   pip install stdb
   pip install geographiclib


Installing development branch from GitHub
-----------------------------------------

.. sourcecode:: bash

   pip install orientpy@git+https://github.com/nfsi-canada/orientpy

Installing from source
----------------------

- Clone the repository:

.. sourcecode:: bash

   git clone https://github.com/nfsi-canada/OrientPy
   cd OrientPy

- Install using pip:

.. sourcecode:: bash

   pip install .
