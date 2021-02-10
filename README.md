
groopsio - Read and write GROOPS file formats in  Python
========================================================

What is groopsio
----------------

groopsio is a free Open Source software package that enables I/O for GROOPS file formats in Python.
It is written as a C/C++ extension for NumPy and requires a few GROOPS  source files.

Installation
------------

To install the current development version of the package, first clone the repository or download the zip archive.
In the root directory of the package (i.e. the directory containing the ``setup.py`` file), run

    GROOPS_SOURCE_DIR=/path/to/groops/source pip install .

to install the package and its dependencies. The environment variable ``GROOPS_SOURCE_DIR`` should point to the source
folder of GROOPS (if ommitted, ``$HOME/groops/source`` is used).
Note that you will need a toolchain capable of building NumPy extensions.
See the [NumPy documentation](https://numpy.org/doc/stable/user/building.html) for details.

License
-------

groopsio a free Open Source software released under the GPLv3 license.


