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
import typing


def loadmat(file_name: str) -> np.ndarray:
    warnings.warn("'loadmat' will be deprecated in favor of 'loadmatrix' in a future release", category=DeprecationWarning)
    return loadmatrix(file_name)


def savemat(file_name: str, M: np.ndarray, mtype: str = 'general', uplo: str = 'upper') -> None:
    warnings.warn("'savemat' will be deprecated in favor of 'savematrix' in a future release", category=DeprecationWarning)
    if uplo.lower() not in ('upper', 'lower'):
        raise ValueError("Matrix triangle must be 'upper' or 'lower'.")

    savematrix(file_name, M, mtype, uplo.lower() == 'lower')


def loadmatrix(file_name: str) -> np.ndarray:
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
    mat : ndarray
        2d ndarray containing the matrix data

    Raises
    ------
    FileNotFoundError
        if file is nonexistent

    Examples
    --------
     >>> import groopsio as gio
     >>> A = gio.loadmatrix('A.dat')
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    return giocpp.loadmat(file_name)


def savematrix(file_name: str, matrix: np.typing.ArrayLike, matrix_type: str = 'general', lower: bool = False) -> None:
    """
    Write Numpy ndarray to GROOPS Matrix file

    Save numeric array in GROOPS matrix file format. Per default
    the array is saved as general matrix. To account for matrix
    properties the keyword arguments mtype and uplo can be used.

    Parameters
    ----------
    file_name : str
        file name
    matrix : array_like(m, n)
        2d ndarray to be written to file
    matrix_type : str
        matrix type {'general', 'symmetric', 'triangular'}, (default: 'general')
    lower : bool
        choose which triangle is stored, ignored when matrix_type='general' (default: the upper triangle is stored)

    Raises
    ------
    FileNotFoundError
        if directory is nonexistent or not writeable

    Examples
    --------
    >>> import numpy as np
    >>> import groopsio as gio
    >>> A = np.eye(10)  # 10x10 identity matrix
    >>>
    >>> gio.savematrix('A.dat', A) # A is saved as general 10x10 matrix
    >>> gio.savematrix('A.dat', A, matrix_type='symmetric') # A is saved as symmetric 10x10 matrix
    >>> gio.savematrix('A.dat', A, matrix_type='triangular', lower=True) # A is saved as lower triangular 10x10 matrix

    """
    matrix = np.asarray(matrix)
    if matrix.ndim == 0:
        warnings.warn('0-dimensional array treated as 1x1 matrix.')
        M = np.atleast_2d(matrix)

    elif matrix.ndim == 1:
        warnings.warn('1-dimensional array treated as column vector.')
        matrix = matrix[:, np.newaxis]

    elif matrix.ndim > 2:
        raise ValueError('ndarray must have at most two dimensions (has {0:d}).'.format(M.ndim))

    if split(file_name)[0] and not isdir(split(file_name)[0]):
        raise FileNotFoundError('Directory ' + split(file_name)[0] + ' does not exist.')

    if matrix_type.lower() not in ('general', 'symmetric', 'triangular'):
        raise ValueError("Matrix type must be 'general', 'symmetric' or 'triangular'.")

    giocpp.savemat(file_name, matrix, matrix_type.lower(), lower)


def loadgridrectangular(file_name: str) -> typing.Tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    """
    Read GROOPS GriddedDataRectangular file format.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    data : (m, n, k) or (m, n) ndarray
        ndarray containing the grid values (if only one value is in the grid file, the last dimension is squeezed)
    lon : (n,) ndarray
        1d ndarray containing the longitude values in radians
    lat : (m,) ndarray
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
    >>> import groopsio.io as gio
    >>> data, lon, lat, a, f = gio.loadgridrectangular('grids/aod1b_RL04.dat')

    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    return_tuple = giocpp.loadgridrectangular(file_name)
    data_count = len(return_tuple) - 4
    if data_count > 1:
        data = np.dstack(return_tuple[0:data_count])
    else:
        data = return_tuple[0]

    return data, return_tuple[-4].flatten(), return_tuple[-3].flatten(), return_tuple[-2], return_tuple[-1]


def loadgrid(file_name: str) -> typing.Tuple[np.ndarray, float, float]:
    """
    Read GROOPS GriddedData file format.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    data : (m, n) ndarray
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
    >>> import groopsio.io as gio
    >>> data, a, f = gio.loadgrid('grids/aod1b_RL04.dat')
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    data, a, f = giocpp.loadgrid(file_name)

    return data, a, f


def savegrid(file_name: str, data: np.typing.ArrayLike, a: float = 6378137.0, f: float = 298.2572221010**-1) -> None:
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


def loadinstrument(file_name: str, concat_arcs: bool = False) -> typing.Tuple[typing.Union[np.ndarray, typing.List[np.ndarray]], int]:
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


def loadinstrumentgnssreceiver(file_name: str) -> typing.Tuple[typing.Dict[str, np.ndarray]]:
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


def loadstarcamera(file_name: str) -> typing.Tuple[np.ndarray, typing.Tuple[np.ndarray]]:
    """
    Read rotation matrices from StarCameraFile.

    Parameters
    ----------
    file_name : str
        file name

    Returns
    -------
    times : (m,) ndarray
        time stamps in MJD
    data : tuple of (3,3) ndarray
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


def saveinstrument(file_name: str, arcs: typing.Union[typing.List[np.ndarray], np.ndarray], epoch_type: typing.Optional[int] = None):
    """
    Save arcs to  GROOPS Instrument file format.

    Parameters
    ----------
    file_name : str
        file name
    arcs : ndarray or list of ndarray
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

    epoch_type = arcs[0].shape[1] - 1 if epoch_type is None else epoch_type

    giocpp.saveinstrument(file_name, [arc for arc in arcs], epoch_type)


def loadsphericalharmonics(file_name: str, retur_sigmas: bool = False) -> typing.Tuple[float, float, np.ndarray, np.ndarray]:
    """
    Read spherical harmonics data from file

    Parameters
    ----------
    file_name : str
        file name
    return_sigmas : bool
        flag whether whether to return formal or calibrated errors (if present)

    Returns
    -------
    anm : array_like(nmax+1, nmax+1)
        Potential coefficients as ndarray. cosine coefficients are stored in the lower triangle, sine coefficients
        above the superdiagonal
    GM : float
        Geocentric gravitational constant
    R : float
        Reference radius
    sigma2anm : array_like(nmax+1, nmax+1)
        variances of potential coefficients in the same structure as anm (if the file does not contain formal erros, a NaN array is returned), only returned if return_sigmas=True.

    Raises
    ------
    FileNotFoundError
        if file is nonexistent
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    return giocpp.loadsphericalharmonics(file_name, return_sigmas)


def savesphericalharmonics(file_name: str, anm: np.typing.ArrayLike, GM: float, R: float, sigma2anm: typing.Optional[np.typing.ArrayLike] = None) -> None:
    """
    Write spherical harmonics data to file

    Parameters
    ----------
    file_name : str
        file name
    anm : array_like(nmax+1, nmax+1)
        Potential coefficients as ndarray. cosine coefficients are stored in the lower triangle, sine coefficients
        above the superdiagonal
    GM : float
        Geocentric gravitational constant
    R : float
        Reference radius
    sigma2anm : array_like(nmax+1, nmax+1)
        Variances of potential coefficients, in the same structure as anm. Default behavior is to not save accuracies (sigma2anm = None).

    Raises
    ------
    FileNotFoundError
        if directory is nonexistent or not writeable

    """
    if split(file_name)[0] and not isdir(split(file_name)[0]):
        raise FileNotFoundError('Directory ' + split(file_name)[0] + ' does not exist.')

    if isinstance(sigma2anm, np.ndarray) and np.all(np.isnan(sigma2anm)):
        sigma2anm = None

    giocpp.savesphericalharmonics(file_name, anm, GM, R, sigma2anm)


def loadgravityfield(file_name):
    warnings.warn("'loadgravityfield' will be deprecated in favor of 'loadsphericalharmonics' in a future release", category=DeprecationWarning)
    anm, GM, R, sigma2anm = giocpp.loadsphericalharmonics(file_name, True)
    return GM, R, anm, sigma2anm


def savegravityfield(file_name, GM, R, anm, sigma2anm=None):
    warnings.warn("'savegravityfield' will be deprecated in favor of 'savesphericalharmonics' in a future release", category=DeprecationWarning)
    savesphericalharmonics(file_name, anm, GM, R, sigma2anm)


def loadtimesplines(file_name: str, time: typing.Union[dt.datetime, float]) -> typing.Tuple[float, float, np.ndarray]:
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
    anm : array_like(nmax+1, nmax+1)
        Potential coefficients as ndarray. cosine coefficients are stored in the lower triangle, sine coefficients
        above the superdiagonal
    GM : float
        Geocentric gravitational constant
    R : float
        Reference radius

    Raises
    ------
    FileNotFoundError
        if file is nonexistent
    """
    if not isfile(file_name):
        raise FileNotFoundError('File ' + file_name + ' does not exist.')

    if isinstance(time, dt.datetime):
        delta = time - dt.datetime(1858, 11, 17)
        time = delta.days + delta.seconds / 86400.0

    GM, R, anm = giocpp.loadtimesplines(file_name, time)

    return anm, GM, R


def loadnormalsinfo(file_name: str, return_full_info: bool = False) -> typing.Tuple[np.ndarray, int, typing.Tuple[str], typing.Optional[np.ndarray], typing.Optional[np.ndarray]]:
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
        raise FileNotFoundError('File ' + splitext(file_name)[0] + '.info.xml' + ' does not exist.')

    lPl, obs_count, names, block_index, used_blocks = giocpp.loadnormalsinfo(file_name)

    if return_full_info:
        return lPl, obs_count, names, block_index.flatten().astype(int), used_blocks.astype(bool)
    else:
        return lPl, obs_count, names


def loadnormals(file_name: str) -> typing.Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
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
        raise FileNotFoundError('File ' + splitext(file_name)[0] + '.info.xml' + ' does not exist.')

    return giocpp.loadnormals(file_name)


def savenormals(file_name: str, N: np.typing.ArrayLike, n: np.typing.ArrayLike, lPl: np.typing.ArrayLike, obs_count: int) -> None:
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
    if (np.asarray(N).ndim != 2) or (np.asarray(N).shape[0] != np.asarray(N).shape[1]):
        raise ValueError('Square normal equation coefficient matrix required.')

    if (np.asarray(n).ndim != 2) or (np.asarray(n).shape[0] != np.asarray(n).shape[0]):
        raise ValueError('Number of parameters in normal equation coefficient matrix and right hand side do not match.')

    if np.asarray(lPl).size != np.asarray(n).shape[1]:
        raise ValueError('Number of right hand sides in observation square sum and right hand side vector do not match.')

    giocpp.savenormals(file_name, N, n, lPl, obs_count)


def loadarclist(file_name: str) -> typing.Tuple[typing.Tuple[int], typing.Tuple[float]]:
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


def loadtimeseries(file_name: str) -> typing.Tuple[np.ndarray, np.ndarray]:
    """
    Read Time Series from matrix/instrument file (based on loadmatrix)

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


def loadfilter(file_name: str) -> typing.Tuple[np.ndarray, np.ndarray, int]:
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


def savefilter(file_name: str, b: np.typing.ArrayLike, a: np.typing.ArrayLike = np.ones(1), start_index: int = 0) -> None:
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

    a = np.asarray(a)
    b = np.asarray(b)

    idx = np.arange(-start_index, -start_index + max(b.size, a.size))

    A = np.zeros((idx.size, 3))
    A[:, 0] = idx
    A[0:b.size, 1] = b
    A[start_index:start_index + a.size, 2] = a

    giocpp.savemat(file_name, A, 0, 0)


def loadpolygon(file_name: str) -> typing.Tuple[np.ndarray]:
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


def savepolygon(file_name: str, polygons: typing.Union[np.typing.ArrayLike, typing.Container[np.typing.ArrayLike]]) -> None:
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


def loadparameternames(file_name: str, encoding: str = 'utf-8', errors: str = 'strict') -> typing.Tuple[str]:
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
