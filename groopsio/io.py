# groopsio - Python wrappers for GROOPS file in-/output.
#
# Copyright (C) 2020 - 2021 GROOPS Developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
"""
Python wrappers for GROOPS file in-/output.
"""

import numpy as np
from os.path import isfile, split, isdir, splitext
import warnings
import datetime as dt
import groopsiobase as giocpp


def loadmat(file_name):
    """
    Read GROOPS Matrix file format.

    Imports GROOPS matrices as numpy array. Note that
    triangular and symmetric matrices are stored fully.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    mat : array_like(m, n)
        2d ndarray containing the matrix data

    Raises
    ------
    FileNotFoundError
        if file is nonexistent

    Examples
    --------
     >>> import groopsio.io as gio
     >>> A = gio.loadmat('A.dat')

    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    return giocpp.loadmat(file_name)


def savemat(file_name, M, mtype='general', uplo='upper'):
    """
    Write Numpy ndarray to GROOPS Matrix file

    Save numeric array in GROOPS matrix file format. Per default
    the array is saved as general matrix. To account for matrix
    properties the keyword arguments mtype and uplo can be used.

    Parameters
    ----------
    file_name : str
        file name
    M : array_like(m, n)
        2d ndarray to be written to file
    mtype : str
        matrix type {'general', 'symmetric', 'triangular'}, (default: 'general')
    uplo : str
        chose which triangle is stored (only applies if mtype is not 'general') {'upper', 'lower'} (default: 'upper')

    Raises
    ------
    FileNotFoundError
        if directory is nonexistent or not writeable

    Examples
    --------
    >>> import numpy as np
    >>> import groopsio.io as gio
    >>> A = np.eye(10)  # 10x10 identity matrix
    >>>
    >>> gio.savemat('A.dat', A) # A is saved as general 10x10 matrix
    >>> gio.savemat('A.dat', A, mtype='symmetric') # A is saved as symmetric 10x10 matrix
    >>> gio.savemat('A.dat', A, mtype='triangular', uplo='lower') # A is saved as lower triangular 10x10 matrix

    """
    if M.ndim == 0:
        warnings.warn('0-dimensional array treated as 1x1 matrix.')
        M = np.atleast_2d(M)

    elif M.ndim == 1:
        warnings.warn('1-dimensional array treated as column vector.')
        M = M[:, np.newaxis]

    elif M.ndim > 2:
        raise ValueError('ndarray must have at most two dimensions (has {0:d}).'.format(M.ndim))

    if split(file_name)[0] and not isdir(split(file_name)[0]):
        raise FileNotFoundError('Directory ' + split(file_name)[0] + ' does not exist.')

    if mtype.lower() not in ('general', 'symmetric', 'triangular'):
        raise ValueError("Matrix type must be 'general', 'symmetric' or 'triangular'.")

    if uplo.lower() not in ('upper', 'lower'):
        raise ValueError("Matrix triangle must be 'upper' or 'lower'.")

    giocpp.savemat(file_name, M, mtype.lower(), uplo.lower())


def loadgridrectangular(file_name):
    """
    Read GROOPS GriddedDataRectangular file format.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    data : list of array_like(m, n)
        2d ndarray containing the grid values
    lon : array_like(n,)
        1d ndarray containing the longitude values in radians
    lat : array_like(m,)
        1d ndarray containing the latitude values in radians
    a : float
        semi-major axis of ellipsoid
    f : float
        flattening of ellipsoid

    Raises
    ------
    FileNotFoundError
        if file is nonexistent

    Examples
    --------
    >>> import numpy as np
    >>> import groopsio.io as gio
    >>> data, a, f = gio.loadgrid('grids/aod1b_RL04.dat')

    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    return_tuple = giocpp.loadgridrectangular(file_name)
    data_count = len(return_tuple) - 4
    return list(return_tuple[0:data_count]), return_tuple[-4].flatten(), return_tuple[-3].flatten(), return_tuple[-2], return_tuple[-1]


def loadgrid(file_name):
    """
    Read GROOPS GriddedData file format.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    data : array_like(m, n)
        2d ndarray containing the grid coordinates and values. Columns 0-3 contain geometry (lon, lat, h, area),
        columns 4-(n-1) contain the corresponding point values
    a : float
        semi-major axis of ellipsoid
    f : float
        flattening of ellipsoid

    Raises
    ------
    FileNotFoundError
        if file is nonexistent

    Examples
    --------
    >>> import numpy as np
    >>> import groopsio.io as gio
    >>> data, a, f = gio.loadgrid('grids/aod1b_RL04.dat')

    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    data, a, f = giocpp.loadgrid(file_name)

    return data, a, f


def savegrid(file_name, data, a=6378137.0, f=298.2572221010**-1):
    """
    Write grid to GROOPS GriddedData file

    Parameters
    ----------
    file_name : str
        file name
    data : array_like(m, n)
        2d ndarray containing the grid coordinates and values. Columns 0-3 contain geometry (lon, lat, h, area),
        columns 4-(n-1) contain the corresponding point values
    a : float
        semi-major axis of ellipsoid
    f : float
        flattening of ellipsoid

    Raises
    ------
    FileNotFoundError
        if directory is nonexistent or not writeable

    Examples
    --------
    >>> import groopsio.io as gio
    >>> G, a, f = gio.loadgrid('grids/aod1b_RL04.dat')
    >>> # manipulate grid
    >>> gio.savegrid('grids/aod1b_RL04_mod.dat', G, a, f)

    """
    if split(file_name)[0] and not isdir(split(file_name)[0]):
        raise FileNotFoundError('Directory ' + split(file_name)[0] + ' does not exist.')

    giocpp.savegrid(file_name, data, a, f)


def loadinstrument(file_name, concat_arcs=False):
    """
    Read GROOPS Instrument file format.

    Instrument data is returned as saved in the file, time stamps are given in MJD.

    Parameters
    ----------
    file_name : str
        file name
    concat_arcs : bool
        flag whether to concatenate all arcs (default: False)

    Returns
    -------
    arcs : tuple of array_like(m, n) or array_like(m, n)
        tuple of 2d ndarrays containing the arc data or a single ndarray if concate_arcs=True
    epoch_Type : int
        enum of instrument type

    Raises
    ------
    FileNotFoundError
        if file is nonexistent

    Examples
    --------
    >>> import numpy as np
    >>> import groopsio.io as gio
    >>> pod2, pod2_type = gio.loadinstrument('satellite/grace2_pod_2008-05.dat')
    >>> pod2, pod2_type = gio.loadinstrument('satellite/grace2_pod_2008-05.dat', concat_arcs=True)

    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    arcs, epoch_type = giocpp.loadinstrument(file_name)

    if concat_arcs:
        arcs = np.hstack(arcs)

    return arcs, epoch_type


def loadinstrumentgnssreceiver(file_name):
    """
    Read GROOPS GnssReceiver (observation or residual) instrument file format.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    arcs : tuple
        A tuple of dictionaries. In each tuple entry, the keys are the combined GROOPS GnssTypes (observation type + satellite type)
        and an additional key for the epochs.
        Keys: epochs , <RINEX signal name><RINEX satellite PRN><(GLONASS) encoded frequency number>, e.g. L1CG10, C1PR09F
        In case of residual files, redundancy and sigma/sigma0 factor for each observation type are added to the dict
        as additional keys '<key>_redundancy' and '<key>_sigmaFactor'. Azimuth and elevation of the residuals are encoded
        as separate 'A1*', 'E1*' (at receiver) and 'A2*', 'E2*' (at transmitter) keys.
        In case multiple entries for one signal exist in one epoch, x characters are appended to the keys of
        the additional observations (e.g. C1CE01, C1CE01x, C1CE01xx, ...).

    Raises
    ------
    FileNotFoundError
        if file is nonexistent


    Examples
    --------
    >>> import numpy as np
    >>> import groopsio.io as gio
    >>> arcs = gio.loadinstrumentgnssreceiver('gnss/gnssReceiver.graz.dat')
    >>> epochs = arcs[0]['epochs']
    >>> obs = arcs[0]['C1CG04']
    """
    if not isfile(file_name):
        raise FileNotFoundError("File {} does not exist.".format(file_name))
    arcs = giocpp.loadinstrumentgnssreceiver(file_name)

    t0 = dt.datetime(1858, 11, 17)
    for arc in arcs:
        is_residual_file = np.any([s.endswith('_redundancy') for s in arc.keys()])

        arc["epochs"] = np.array([t0 + dt.timedelta(days=tk) for tk in arc["epochs"]])
        to_delete = []
        # Remove 0 values in observations
        for obs_name, values in arc.items():
            if obs_name.endswith("_redundancy") or obs_name.endswith("_sigmaFactor") or obs_name[0] == "D" or '*' in obs_name[0:3] or obs_name == "epochs":
                continue
            if np.count_nonzero(values[~np.isnan(values)]) == 0:
                to_delete.extend([obs_name, obs_name + '_redundancy', obs_name + '_sigmaFactor'] if is_residual_file else [obs_name])

        for obs_name in to_delete:
            arc.pop(obs_name)

    return arcs


def loadstarcamera(file_name):
    """
    Read rotation matrices from StarCameraFile.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    times : array_like(m,)
        time stamps in MJD
    data : tuple of array_like(3,3)
        rotation matrices for each epoch

    Raises
    ------
    FileNotFoundError
        if file is nonexistent
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    times, data = giocpp.loadstarcamera(file_name)

    return times.flatten(), data


def saveinstrument(file_name, arcs, epoch_type=None):
    """
    Save arcs to  GROOPS Instrument file format.

    Parameters
    ----------
    file_name : str
        file name
    arcs : list of array_like(m, n) or array_like(m, n)
        arc-wise data as ndarray, or single ndarray
    epoch_type : int
        enum of epoch type (Default: MISCVALUES)

    Raises
    ------
    FileNotFoundError
        if directory is nonexistent or not writeable

    Examples
    --------
    >>> import numpy as np
    >>> import groopsio.io as gio
    >>> pod2, pod2_type = gio.loadinstrument('satellite/grace2_pod_2008-05.dat')
    >>> gio.saveinstrument('tmp/grace2_pod_2008-05_arcs_1-5-17.dat', pod2, pod2_type)

    """
    if split(file_name)[0] and not isdir(split(file_name)[0]):
        raise FileNotFoundError('Directory ' + split(file_name)[0] + ' does not exist.')

    if type(arcs) is not list:
        arcs = [arcs]

    epoch_type = arcs[0].shape[1]-1 if epoch_type is None else epoch_type

    giocpp.saveinstrument(file_name, [arc for arc in arcs], epoch_type)


def loadgravityfield(file_name):
    """
    Read SphericalHarmonics from gfc-file

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    GM : float
        Geocentric gravitational constant
    R : float
        Reference radius
    anm : array_like(nmax+1, nmax+1)
        Potential coefficients as ndarray. cosine coefficients are stored in the lower triangle, sine coefficients
        above the superdiagonal
    sigma2anm : array_like(nmax+1, nmax+1)
        Variances of potential coefficients, in the same structure as anm. If the gfc file does not provide accuracies,
        a NAN array is returned.

    Raises
    ------
    FileNotFoundError
        if file is nonexistent
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    GM, R, anm, sigma2anm = giocpp.loadgravityfield(file_name)

    return GM, R, anm, sigma2anm


def savegravityfield(file_name, GM, R, anm, sigma2anm=None):
    """
    Write GravityField instance to gfc-file

    Parameters
    ----------
    file_name : str
        file name
    GM : float
        Geocentric gravitational constant
    R : float
        Reference radius
    anm : array_like(nmax+1, nmax+1)
        Potential coefficients as ndarray. cosine coefficients are stored in the lower triangle, sine coefficients
        above the superdiagonal
    sigma2anm : array_like(nmax+1, nmax+1)
        Variances of potential coefficients, in the same structure as anm. Default behavior is to not save accuracies
        (sigma2anm = None).

    Raises
    ------
    FileNotFoundError
        if directory is nonexistent or not writeable

    """
    if split(file_name)[0] and not isdir(split(file_name)[0]):
        raise FileNotFoundError('Directory ' + split(file_name)[0] + ' does not exist.')

    has_sigmas = sigma2anm is not None

    giocpp.savegravityfield(file_name, GM, R, anm, has_sigmas, sigma2anm if has_sigmas else None)


def loadtimesplines(file_name, time):
    """
    Read potential coefficients from TimeSplines file


    Parameters
    ----------
    file_name : str
        file name
    time : float or datetime.datetime
        evaluation time of TimeSplines file as MJD (float) or datetime object

    Returns
    -------
    GM : float
        Geocentric gravitational constant
    R : float
        Reference radius
    anm : array_like(nmax+1, nmax+1)
        Potential coefficients as ndarray. cosine coefficients are stored in the lower triangle, sine coefficients
        above the superdiagonal

    Raises
    ------
    FileNotFoundError
        if file is nonexistent
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    if isinstance(time, dt.datetime):
        delta = time-dt.datetime(1858, 11, 17)
        time = delta.days + delta.seconds/86400.0

    GM, R, anm = giocpp.loadtimesplines(file_name, time)

    return GM, R, anm


def loadnormalsinfo(file_name, return_full_info=False):
    """
    Read metadata of normal equation file.

    Parameters
    ----------
    file_name : str
        file name of normal equations
    return_full_info : bool
        if true, return lPl, observation count, parameter names, block index and used blocks, else (default)
        return only lPl, observation count and parameter names

    Returns
    -------
    lPl : array_like(rhs_count,)
        square sum of observations for each right hand side
    obs_count : int
        observation count
    names : tuple of str
        string representation of parameter names
    block_index : array_like(block_count+1,)
        beginning/end of normal equation blocks. Only returned if return_full_info is true.
    used_blocks : array_like(block_count, block_count)
        boolean array representing the sparsity structure of the normal equations. Only returned if return_full_info
        is true.

    Raises
    ------
    FileNotFoundError
        if file is nonexistent
    """
    if not isfile(splitext(file_name)[0] + '.info.xml'):
        raise FileNotFoundError('File ' + splitext(file_name)[0] +'.info.xml' + ' does not exist.')

    lPl, obs_count, names, block_index, used_blocks = giocpp.loadnormalsinfo(file_name)

    if return_full_info:
        return lPl, obs_count, names, block_index.flatten().astype(int), used_blocks.astype(bool)
    else:
        return lPl, obs_count, names


def loadnormals(file_name):
    """
    Read normal equations from file file.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    N : array_like(parameter_count, parameter_count)
        coefficient matrix of normal equations
    n : array_like(parameter_count, ŕhs_count)
        right hand side(s)
    lPl : array_like(rhs_count,)
        square sum of observations for each right hand side
    obs_count : int
        observation count

    Raises
    ------
    FileNotFoundError
        if file is nonexistent
    """
    if not isfile(splitext(file_name)[0] + '.info.xml'):
        raise FileNotFoundError('File ' + splitext(file_name)[0] +'.info.xml' + ' does not exist.')

    return giocpp.loadnormals(file_name)


def savenormals(file_name, N, n, lPl, obs_count):
    """
    Read normal equations from file file.

    Parameters
    ----------
    file_name : str
        file name
    N : array_like(m, m)
        normal equation matrix
    n : array_like(m, k)
        right hand side
    lPl : array_like(k,)
        square sum of observations
    obs_count : int
        number of observations

    Raises
    ------
    ValueError
        if dimensions of passed arguments do not match

    """
    if (N.ndim != 2) or (N.shape[0] != N.shape[1]):
        raise ValueError('Square normal equation coefficient matrix required.')

    if (n.ndim != 2) or (n.shape[0] != N.shape[0]):
        raise ValueError('Number of parameters in normal equation coefficient matrix and right hand side do not match.')

    if lPl.size != n.shape[1]:
        raise ValueError('Number of right hand sides in observation square sum and right hand side vector do not match.')

    return giocpp.savenormals(file_name, N, n, lPl, obs_count)


def loadarclist(file_name):
    """
    Read GROOPS arcList file.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    arc_intervals : tuple of int
        interval bounds in arc indices
    time_intervals : tuple of float
        interval bounds in MJD

    Raises
    ------
    FileNotFoundError
        if file does not exist
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    return giocpp.loadarclist(file_name)


def loadtimeseries(file_name):
    """
    Read Time Series from matrix/instrument file (based on loadmat)

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    t : array_like(m,)
        time stamps in MJD
    X : array_like(m, n)
        data values associated with each time stamp
    Raises
    ------
    FileNotFoundError
        if file not found
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    ts = loadmat(file_name)

    return ts[:, 0], ts[:, 1::]


def loadfilter(file_name):
    """
    Read digital filter from file.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    b : array_like(p,)
        MA coefficients
    a : array_like(q,)
        AR coefficients
    start_index : int
        positive integer which determines the time shift (non-causality) of the MA coefficients

    Raises
    ------
    FileNotFoundError
        if file does not exist
    ValueError
        if a negative index of an AR coefficient is found (AR part must be causal)
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    A = giocpp.loadmat(file_name)
    idx_bk = A[A[:, 1] != 0, 0].astype(int)
    idx_ak = A[A[:, 2] != 0, 0].astype(int)

    if np.any(idx_ak < 0):
        raise ValueError('Negative indices for AR coefficients not allowed (causal-filter).')

    start_index = max(-np.min(idx_bk), 0)

    b = np.zeros(np.max(idx_bk) + start_index + 1)
    b[idx_bk + start_index] = A[A[:, 1] != 0, 1]
    a = np.zeros(np.max(idx_ak) + 1 if idx_ak.size > 0 else 1)
    a[idx_ak] = A[A[:, 2] != 0, 2]
    if a[0] != 1.0:
        a[0] = 1.0
        warnings.warn('a0 coefficient set to one.')

    return b, a, start_index


def savefilter(file_name, b, a=np.ones(1), start_index=0):
    """
    Save filter coefficients to file.

    Parameters
    ----------
    file_name : str
        file name
    b : array_like(p,)
        MA coefficients
    a = array_like(q,)
        AR coefficients (Default: [1], pure MA filter)
    start_index : int
        positive integer which determines the time shift (non-causality) of the MA coefficients (Default: 0)

    Raises
    ------
    ValueError
        if a negative start_index is passed
    FileNotFoundError
        if directory is not writeable
    """
    if start_index < 0:
        raise ValueError('start_index must be positive')

    if split(file_name)[0] and not isdir(split(file_name)[0]):
        raise FileNotFoundError('Directory ' + split(file_name)[0] + ' does not exist.')

    idx = np.arange(-start_index, -start_index + max(b.size, a.size))
    print(idx)

    A = np.zeros((idx.size, 3))
    A[:, 0] = idx
    A[0:b.size, 1] = b
    A[start_index:start_index + a.size, 2] = a

    giocpp.savemat(file_name, A, 0, 0)


def loadpolygon(file_name):
    """
    Read  a polygon list from file.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    polygons : tuple of array_likes(p, 2)
        tuple of 2-d arrays representing the vertices (longitude, latitude) of each polygon in radians.

    Raises
    ------
    FileNotFoundError
        if file does not exist
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    return giocpp.loadpolygon(file_name)


def savepolygon(file_name, polygons):
    """
    Save a polygon list to file.

    Parameters
    ----------
    file_name : str
        file name
    polygons : array_like(p, 2) or tuple of array_likes(p, 2)
        tuple of 2-d arrays representing the vertices (longitude, latitude) of each polygon in radians.

    Raises
    ------
    FileNotFoundError
        if directory is not writeable
    """
    if split(file_name)[0] and not isdir(split(file_name)[0]):
        raise FileNotFoundError('Directory ' + split(file_name)[0] + ' does not exist.')

    if type(polygons) is list:
        polygons = tuple(polygons)

    if type(polygons) is not tuple:
        polygons = (polygons,)

    return giocpp.savepolygon(file_name, polygons)


def loadparameternames(file_name, encoding='utf-8', errors='strict'):
    """
    Read a parameter name list from file.

    Parameters
    ----------
    file_name : str
        file name
    encoding : str
        encoding used to decode the bytes (see bytes.decode)
    error : str
        error level for decoding (see bytes.decode)

    Returns
    -------
    parameter_names : tuple of str
        tuple of strings

    Raises
    ------
    FileNotFoundError
        if file does not exist
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    return tuple(name.decode(encoding, errors) for name in giocpp.loadparameternames(file_name))
