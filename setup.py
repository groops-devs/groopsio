import setuptools
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
import numpy as np
import os


def configuration(parent_package='', top_path=None):

    groops_dir = os.getenv('GROOPS_SOURCE_DIR')
    if groops_dir is None:
        groops_dir = os.path.join(os.path.expanduser('~'), 'groops', 'source')

    include_dirs = [groops_dir]
    library_dirs = []
    libraries = ['expat', 'z', 'gfortran', 'stdc++fs']

    source_files = []
    with open('sources.list', 'r') as f:
        for line in f:
            if len(line) > 0 and not line.startswith('#'):
                source_files.append(os.path.join(groops_dir, line.strip()))

    lapack_opts = np.__config__.lapack_opt_info
    try:
        include_dirs.extend(lapack_opts['include_dirs'])
    except KeyError:
        pass
    try:
        library_dirs.extend(lapack_opts['library_dirs'])
    except KeyError:
        pass
    try:
        libraries.extend(lapack_opts['libraries'])
    except KeyError:
        pass
    blas_opts = np.__config__.blas_opt_info
    try:
        include_dirs.extend(blas_opts['include_dirs'])
    except KeyError:
        pass
    try:
        library_dirs.extend(blas_opts['library_dirs'])
    except KeyError:
        pass
    try:
        libraries.extend(blas_opts['libraries'])
    except KeyError:
        pass

    config = Configuration()
    config.add_installed_library('groopsdeps', source_files, 'groopsio/lib',
                                 build_info={'include_dirs': include_dirs, 'libraries': libraries,
                                 'library_dirs': library_dirs})

    libraries.append('groopsdeps')
    library_dirs.append('groopsio/lib')
    config.add_extension('groopsiobase',  ["src/groopsio.cpp"], language='c++',
                        include_dirs=include_dirs, libraries=libraries, library_dirs=library_dirs)

    return config


setup(
    name='groopsio',
    version='0.1',
    author='Andreas Kvas',
    description='A python package to enable I/O for GROOPS files',
    packages=['groopsio'],
    configuration=configuration
)
