
groopsio - Read and write GROOPS file formats in  Python
========================================================

What is groopsio
----------------

groopsio is a free Open Source software package for reading and writing GROOPS file formats in Python.
It is written as a C/C++ extension to NumPy and requires a few GROOPS  source files.

Installation
------------

To install the current development version of the package, first clone the repository or download the zip archive.
In the root directory of the package (i.e. the directory containing the ``setup.py`` file), run

    GROOPS_SOURCE_DIR=/path/to/groops/source pip install .

to install the package and its dependencies. The environment variable ``GROOPS_SOURCE_DIR`` should point to the source
folder of GROOPS (if ommitted, ``$HOME/groops/source`` or the Windows equivalent is used).

Note that you will need a toolchain capable of building GROOPS. See the GROOPS installation guide for
[Windows](https://github.com/groops-devs/groops/blob/main/INSTALL.md#microsoft-windows)
or [Linux](https://github.com/groops-devs/groops/blob/main/INSTALL.md#linux) for details (no optional steps are necessary for groopsio).

License
-------

groopsio a free Open Source software released under the GPLv3 license.
