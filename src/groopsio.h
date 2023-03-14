/*
 * Python module for GROOPS file I/O.
 *
 * Copyright (C) 2020 - 2021 GROOPS Developers
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#ifndef __GROOPSIO__
#define __GROOPSIO__

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include "numpy/ndarraytypes.h"
#include "numpy/npy_3kcompat.h"

#include "base/matrix.h"
#include "inputOutput/file.h"
#include "parallel/parallel.h"

#include "files/filePlatform.h"
#include "files/fileAdmittance.h"
#include "files/fileArcList.h"
#include "files/fileDoodsonEarthOrientationParameter.h"
#include "files/fileDoodsonHarmonic.h"
#include "files/fileEarthOrientationParameter.h"
#include "files/fileEarthTide.h"
#include "files/fileEphemerides.h"
#include "files/fileGnssAntennaDefinition.h"
#include "files/fileGnssReceiverDefinition.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileGriddedData.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"
#include "files/fileMeanPolarMotion.h"
#include "files/fileNormalEquation.h"
#include "files/fileOceanPoleTide.h"
#include "files/fileParameterName.h"
#include "files/filePolygon.h"
#include "files/fileSatelliteModel.h"
#include "files/fileSphericalHarmonics.h"
#include "files/fileStringTable.h"
#include "files/fileTideGeneratingPotential.h"
#include "files/fileTimeSplinesGravityfield.h"
#include "files/fileVariationalEquation.h"

static PyObject *groopsError;

/*
 * Create 2D Python numeric array of specific size (FORTRAN order).
 */
PyObject* createDoubleArray(UInt rows, UInt columns)
{
  npy_intp size[2];
  size[0] = rows;
  size[1] = columns;

  PyObject *pyObj = PyArray_ZEROS(2, size, NPY_DOUBLE, 1/*fortran order*/);

  return pyObj;
}

/*
 * GROOPS Matrix to Python Array
 */
PyObject* fromMatrix(const Matrix& M)
{
  if(M.getType() == Matrix::SYMMETRIC)
    fillSymmetric(M);
  else if(M.getType() == Matrix::TRIANGULAR)
    zeroUnusedTriangle(M);

  PyObject *pyObj = createDoubleArray(M.rows(), M.columns());
  memcpy(PyArray_DATA((PyArrayObject*)(pyObj)), M.field(), M.size()*sizeof(Double));

  return pyObj;
}

/*
 * Python Array to GROOPS Matrix
 */
Matrix fromPyObject(PyObject *pyObj, Matrix::Type t = Matrix::GENERAL,
                    Matrix::Uplo uplo = Matrix::UPPER)
{
  PyObject *pyArr = PyArray_FROM_OTF(pyObj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
  if(!pyArr)
    throw(Exception("Unable to convert input to float array."));
  Int ndims = PyArray_NDIM((PyArrayObject*)pyArr);
  if(ndims<1 || ndims>2)
    throw(Exception("Expected one or two dimensional Numpy Array (got "+ndims%"%i"s+" dimensions)."));
  npy_intp* size = PyArray_DIMS((PyArrayObject*)pyArr);
  UInt rows = size[0];
  UInt columns = ndims == 1 ? 1 : size[1];

  Matrix M;
  if(t==Matrix::GENERAL)
    M = Matrix(rows, columns);
  else
  {
    if(rows != columns)
      throw(Exception("Matrix must be square."));
    M = Matrix(rows, t, uplo);
  }

  if(size[0]*size[1] != 0)
    memcpy(M.field(), PyArray_DATA((PyArrayObject*)pyArr), M.size()*sizeof(Double));

  Py_DecRef(pyArr);
  return M;
}

/*
 * Read GROOPS matrix from file
 */
static PyObject* loadmatrix(PyObject* /*self*/, PyObject *args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments."));

    Matrix M;
    readFileMatrix(FileName(std::string(s)), M);
    PyObject* pyObj = fromMatrix(M);

    return pyObj;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Save numeric array to GROOPS matrix file format
 */
static PyObject* savematrix(PyObject* /*self*/, PyObject *args)
{
  try
  {
    const char *s;
    const char *type;
    PyObject *pyObj, *isLower;
    if(!PyArg_ParseTuple(args, "sOsO", &s, &pyObj, &type, &isLower))
      throw(Exception("Unable to parse arguments."));

    Matrix::Type type_rep = Matrix::GENERAL;
    if(std::string(type) == std::string("symmetric"))
      type_rep = Matrix::SYMMETRIC;
    else if(std::string(type) == std::string("triangular"))
      type_rep = Matrix::TRIANGULAR;

    Matrix::Uplo uplo_rep = PyObject_IsTrue(isLower) ? Matrix::LOWER : Matrix::UPPER;

    Matrix M = fromPyObject(pyObj, type_rep, uplo_rep);
    writeFileMatrix(FileName(std::string(s)), M);

    Py_RETURN_NONE;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Read GROOPS grid rectangular from file
 */
static PyObject* loadgridrectangular(PyObject* /*self*/, PyObject *args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments."));

    GriddedDataRectangular G;
    readFileGriddedData(FileName(std::string(s)), G);

    const UInt dataCount = G.values.size();

    PyObject *return_tuple = PyTuple_New(dataCount + 4); // data, lon, lat, a, f

    for(UInt k = 0; k < dataCount; k++)
      PyTuple_SetItem(return_tuple, k, fromMatrix(G.values.at(k)));

    Matrix lons(G.longitudes.size(), 1);
    for(UInt k = 0; k < G.longitudes.size(); k++)
      lons(k, 0) = G.longitudes.at(k);

    Matrix lats(G.latitudes.size(), 1);
    for(UInt k = 0; k < G.latitudes.size(); k++)
      lats(k, 0) = G.latitudes.at(k);

    PyTuple_SetItem(return_tuple, dataCount + 0, fromMatrix(lons));
    PyTuple_SetItem(return_tuple, dataCount + 1, fromMatrix(lats));
    PyTuple_SetItem(return_tuple, dataCount + 2, PyFloat_FromDouble(G.ellipsoid.a()));
    PyTuple_SetItem(return_tuple, dataCount + 3, PyFloat_FromDouble(G.ellipsoid.f()));

    return return_tuple;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Read GROOPS grid from file
 */
static PyObject* loadgrid(PyObject* /*self*/, PyObject *args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments."));

    GriddedData G;
    readFileGriddedData(FileName(std::string(s)), G);

    PyObject* data = createDoubleArray(G.points.size(), 4 + G.values.size()); // lon, lat, h, area, values
    {
      Angle L, B;
      Double h;

      for(UInt k=0; k<G.points.size(); k++)
      {
        G.ellipsoid(G.points.at(k), L, B, h);
        *(static_cast<Double*>(PyArray_GETPTR2((PyArrayObject*)data, k, 0))) = static_cast<Double>(L);
        *(static_cast<Double*>(PyArray_GETPTR2((PyArrayObject*)data, k, 1))) = static_cast<Double>(B);
        *(static_cast<Double*>(PyArray_GETPTR2((PyArrayObject*)data, k, 2))) = h;
        if(G.areas.size() == G.points.size()) *(static_cast<Double*>(PyArray_GETPTR2((PyArrayObject*)data, k, 3))) = G.areas.at(k);
        for(UInt l = 0; l<G.values.size(); l++)
          *(static_cast<Double*>(PyArray_GETPTR2((PyArrayObject*)data, k, 4+l))) = G.values.at(l).at(k);
      }
    }

    PyObject *return_tuple = PyTuple_New(3); // data, a, f
    PyTuple_SetItem(return_tuple, 0, data);
    PyTuple_SetItem(return_tuple, 1, PyFloat_FromDouble(G.ellipsoid.a()));
    PyTuple_SetItem(return_tuple, 2, PyFloat_FromDouble(G.ellipsoid.f()));

    return return_tuple;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Save Python objects to GROOPS grid file format
 */
static PyObject* savegrid(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    Double a = 0.0;
    Double f = 0.0;
    PyObject *data_array;
    if(!PyArg_ParseTuple(args, "sOdd", &s, &data_array, &a, &f))
      throw(Exception("Unable to parse arguments."));

    Ellipsoid ell(a, f != 0.0 ? 1/f : 0.0);
    Matrix data = fromPyObject(data_array);
    const UInt pointCount = data.rows();
    const UInt valueCount = data.columns() - 4;

    std::vector<Vector3d> _point(pointCount);
    std::vector<Double> _area(pointCount);
    std::vector< std::vector<Double> > _value(valueCount, std::vector<Double>(pointCount));

    for(UInt k = 0; k<pointCount; k++)
    {
      _point.at(k) = ell(Angle(data(k, 0)), Angle(data(k, 1)), Angle(data(k, 2)));
      _area.at(k) = data(k, 3);
      for(UInt l = 0; l<valueCount; l++)
        _value.at(l).at(k) = data(k, 4+l);
    }

    writeFileGriddedData(FileName(std::string(s)), GriddedData(ell, _point, _area, _value));

    Py_RETURN_NONE;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Load spherical harmonic coefficients from GROOPS file
 */
static PyObject* loadsphericalharmonics(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    PyObject *returnErrors;
    if(!PyArg_ParseTuple(args, "sO", &s, &returnErrors))
      throw(Exception("Unable to parse arguments."));

    SphericalHarmonics coeffs;
    readFileSphericalHarmonics(FileName(std::string(s)), coeffs);

    PyObject *return_tuple = PyTuple_New((PyObject_IsTrue(returnErrors) ? 4 : 3));

    Matrix anm = coeffs.cnm();
    anm.setType(Matrix::GENERAL);
    axpy(1.0, coeffs.snm().slice(1, 1, coeffs.maxDegree(), coeffs.maxDegree()).trans(), anm.slice(0, 1, coeffs.maxDegree(), coeffs.maxDegree()));
    PyTuple_SetItem(return_tuple, 0, fromMatrix(anm));
    PyTuple_SetItem(return_tuple, 1, PyFloat_FromDouble(coeffs.GM()));
    PyTuple_SetItem(return_tuple, 2, PyFloat_FromDouble(coeffs.R()));

    if(PyObject_IsTrue(returnErrors))
    {
      Matrix sigma2anm(coeffs.maxDegree()+1, coeffs.maxDegree()+1, NAN_EXPR);
      if(coeffs.sigma2cnm().size() && ((quadsum(coeffs.sigma2cnm())+quadsum(coeffs.sigma2snm())) != 0))
      {
        sigma2anm = coeffs.sigma2cnm();
        sigma2anm.setType(Matrix::GENERAL);
        axpy(1.0, coeffs.sigma2snm().slice(1, 1, coeffs.maxDegree(), coeffs.maxDegree()).trans(), sigma2anm.slice(0, 1, coeffs.maxDegree(), coeffs.maxDegree()));
      }
      PyTuple_SetItem(return_tuple, 3, fromMatrix(sigma2anm));
    }

    return return_tuple;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Save Python objects to a gfc-file
 */
static PyObject* savesphericalharmonics(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    Double GM = 0.0;
    Double R = 0.0;
    PyObject *anm, *sigma2anm;
    if(!PyArg_ParseTuple(args, "sOddO", &s, &anm, &GM, &R, &sigma2anm))
      throw(Exception("Unable to parse arguments."));

    Bool hasSigmas = sigma2anm != Py_None;

    Matrix _cnm = fromPyObject(anm);
    const UInt maxDegree = _cnm.rows()-1;
    Matrix _snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    copy(_cnm.slice(0, 1, maxDegree, maxDegree).trans(), _snm.slice(1, 1, maxDegree, maxDegree));
    zeroUnusedTriangle(_snm);
    _cnm.setType(Matrix::TRIANGULAR, Matrix::LOWER);
    zeroUnusedTriangle(_cnm);

    SphericalHarmonics harm;
    if(hasSigmas)
    {
      Matrix _sigma2cnm = fromPyObject(sigma2anm);
      Matrix _sigma2snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      copy(_sigma2cnm.slice(0, 1, maxDegree, maxDegree).trans(), _sigma2snm.slice(1, 1, maxDegree, maxDegree));
      zeroUnusedTriangle(_snm);
      _sigma2cnm.setType(Matrix::TRIANGULAR, Matrix::LOWER);
      zeroUnusedTriangle(_sigma2cnm);

      harm = SphericalHarmonics(GM, R, _cnm, _snm, _sigma2cnm, _sigma2snm);
    }
    else
      harm = SphericalHarmonics(GM, R, _cnm, _snm);

    writeFileSphericalHarmonics(FileName(std::string(s)), harm);

    Py_RETURN_NONE;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Load spherical harmonic coefficients from a TimeSplines file
 */
static PyObject* loadtimesplines(PyObject* /*self*/, PyObject* args)
{
  try
  {
    Double mjd = 0.0;
    const char *s;
    if(!PyArg_ParseTuple(args, "sd", &s, &mjd))
      throw(Exception("Unable to parse arguments."));
    std::string fname(s);

    InFileTimeSplinesGravityfield timeSplinesFile;
    timeSplinesFile.open(FileName(fname));

    Time t = mjd2time(mjd);
    SphericalHarmonics coeffs = timeSplinesFile.sphericalHarmonics(t);
    timeSplinesFile.close();

    Matrix anm = coeffs.cnm();
    anm.setType(Matrix::GENERAL);
    axpy(1.0, coeffs.snm().slice(1, 1, coeffs.maxDegree(), coeffs.maxDegree()).trans(), anm.slice(0, 1, coeffs.maxDegree(), coeffs.maxDegree()));

    PyObject *return_tuple = PyTuple_New(3);
    PyTuple_SetItem(return_tuple, 0, PyFloat_FromDouble(coeffs.GM()));
    PyTuple_SetItem(return_tuple, 1, PyFloat_FromDouble(coeffs.R()));
    PyTuple_SetItem(return_tuple, 2, fromMatrix(anm));

    return return_tuple;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Load arcs from instrument files
 */
static PyObject* loadinstrument(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments."));
    std::string fname(s);

    InstrumentFile fileInstrument;
    fileInstrument.open(FileName(fname));

    PyObject* return_values = PyTuple_New(2); // (time data), EpochType

    PyObject* arc_list = PyTuple_New(fileInstrument.arcCount());
    Epoch::Type type = Epoch::EMPTY;

    for(UInt arcNo = 0; arcNo<fileInstrument.arcCount(); arcNo++)
    {
      Arc arc = fileInstrument.readArc(arcNo);

      Matrix M;
      if(arc.size())
      {
        M = arc.matrix();
        type = arc.at(0).getType();
      }

      PyTuple_SetItem(arc_list, arcNo, fromMatrix(M));
    }

    PyTuple_SetItem(return_values, 0, arc_list);
    PyTuple_SetItem(return_values, 1, Py_BuildValue("i", static_cast<Int>(type)));

    return return_values;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Load arcs from GnssReceiver instrument file
 */
static PyObject* loadinstrumentgnssreceiver(PyObject *, PyObject* args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments"));
    std::string fname(s);
    InstrumentFile file_receiver(fname);

    PyObject* return_tuple = PyTuple_New(file_receiver.arcCount());
    const std::string epoch_str = "epochs";
    const std::vector<std::string> repetiton_obs_name_additions = {"_redundancy","_sigmaFactor"};
    for(UInt arcNo = 0; arcNo < file_receiver.arcCount(); arcNo++)
    {
      GnssReceiverArc arc = file_receiver.readArc(arcNo);
      PyObject *arc_dict = PyDict_New();
      if(!arc_dict)
        throw(Exception("could not create dictionary"));

      MiscValuesArc arc_new;
      auto firstEpoch = arc.front().time;
      auto lastEpoch = arc.back().time;
      Vector deltat(arc.size());
      for(UInt i = 0; i < arc.size() - 1; i++)
      {
        deltat(i) = (arc.at(i+1).time.mjd() - arc.at(i).time.mjd()) * 24 * 3600;
      }
      const Double dt = median(deltat);
      UInt current_epoch = 0;
      const UInt samples = static_cast<UInt>(std::round(((lastEpoch.mjd() - firstEpoch.mjd()) * 24 * 3600) / dt)) + 1;

      // Create an array of zeros and inser them into the dict with the key ["epochs"]
      // works a intendet
      npy_intp dims[0];
      dims[0] = samples;
      PyObject *times =  PyArray_ZEROS(1, dims, NPY_DOUBLE, 1);
      if(!times)
        throw(Exception("could not create array"));


      PyDict_SetItemString(arc_dict, epoch_str.c_str() , times);
      for(auto& epoch : arc)
      {
        // GET current epoch and the idx in the 1d array and change the value in the dict
        // is effectively arc_dict["epochs"][idx] = VALUE
        Double mjd = epoch.time.mjd();
        PyArrayObject *m = (PyArrayObject*)PyDict_GetItemString(arc_dict, epoch_str.c_str());
        const Int idx = static_cast<Int>(std::round((mjd - firstEpoch.mjd()) * 24 * 3600 / dt));
        *(static_cast<Double*>(PyArray_GETPTR1(m, idx))) = mjd;
        UInt idObs = 0;
        // go over all satellites in this epoch
        for(GnssType typeSat : epoch.satellite)
        {
          // searach int he residualfile list of obstype when the satellitesystem starts
          // example ***E10 would be distance = 8 for:
          // A1*G** 0
          // E1*G** 1
          // A2*G** 2
          // E2*G** 3
          // I**G** 4
          // I**G** 5
          // I**G** 6
          // C1CG** 7
          // A1*E** 8
          UInt idType = std::distance(epoch.obsType.begin(), std::find(epoch.obsType.begin(), epoch.obsType.end(), typeSat));

          // vector with types that have been used already. Is used to check if same named values as used for redundancy and
          // sigmafactor
          std::vector<std::string> used_types;

          // id obs is initialized for every epoch with 0 -> is the index for the value
          // idType is initialized for every satellite -> is the index for the obsType
          // while current idType is still valid for the current satellite and idObs within the vector of values
          while((idType<epoch.obsType.size()) && (idObs<epoch.observation.size()) && (epoch.obsType.at(idType) == typeSat))
          {
             // create obs and satellite combined gnsstype and iterate idType to next obsType e.g. L1CG10
             GnssType obs = epoch.obsType.at(idType++);
             GnssType full = obs + typeSat;
             std::string str_type = full.str();
             std::string str_type_tmp = str_type;
             // get current value for the fulltype and iteratoe idobs to next observation value
             Double value = epoch.observation.at(idObs++);

             // If the current fullType is within the used types append something to the name string
             // If L1CG10 already in usedTypes will append _redundancy to it leading to -> L1CG10_redundancy
             // if L1CG10_redundancy already in usedTypes will append _sigmaFactor to the name string .> L1CG10_sigmaFactor
             // if L1CG10_sigmaFactor in usedTypes will append "X" to the name until no longer in used types
             UInt repetition = 0;
             while(std::find(used_types.begin(), used_types.end(), str_type) != used_types.end())
             {
               if(repetition < repetiton_obs_name_additions.size())
               {
                 str_type = str_type_tmp + repetiton_obs_name_additions.at(repetition);
               }
               else
               {
                 str_type += "x";
               }
               repetition++;
             }


             used_types.push_back(str_type);
             // if the current full type string is not in the created dict add it to the dict with an array of zeros
             // with the same length as epoch
             if (!PyDict_GetItemString(arc_dict, str_type.c_str()))
             {
                // create and add zero array for the given full type
               PyObject* zVals = PyArray_ZEROS(1, dims, NPY_DOUBLE, 1);
               Int errorCode = PyDict_SetItemString(arc_dict, str_type.c_str(), zVals);
               if( errorCode == -1)
                 std::cout << "Error while emplacing: " << str_type.c_str() << " epoch idx: "  << idx << std::endl;

               PyArrayObject* val_ar_temp = (PyArrayObject*)PyDict_GetItemString(arc_dict, str_type.c_str());
               if(val_ar_temp == NULL)
                 std::cout << "Could not find: " << str_type << " epoch idx: "  << idx << std::endl;

               for(UInt i = 0 ; i < samples; ++i)
                 *(static_cast<Double*>(PyArray_GETPTR1(val_ar_temp, i)))  = NAN_EXPR;

            }
            PyArrayObject *arPointer = (PyArrayObject*)PyDict_GetItemString(arc_dict, str_type.c_str());
            *(static_cast<Double*>(PyArray_GETPTR1(arPointer, idx))) = value;
          }
        }
      }
      PyTuple_SetItem(return_tuple, arcNo, arc_dict);
    }
    return return_tuple;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Load arcs from a StarCamera instrument file
 */
static PyObject* loadstarcamera(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments."));
    std::string fname(s);

    StarCameraArc starCameraArc = InstrumentFile::read(FileName(fname));

    PyObject* return_values = PyTuple_New(2); // time, data

    PyObject* data = PyTuple_New(starCameraArc.size());

    Vector mjd(starCameraArc.size());
    for(UInt k = 0; k<starCameraArc.size(); k++)
    {
      Time epoch = starCameraArc.at(k).time;
      Rotary3d rot = starCameraArc.at(k).rotary;

      PyTuple_SetItem(data, k, fromMatrix(rot.matrix()));
      mjd(k) = epoch.mjd();
    }

    PyTuple_SetItem(return_values, 0, fromMatrix(mjd));
    PyTuple_SetItem(return_values, 1, data);

    return return_values;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Save a list of numpy arrays as instrument file
 */
static PyObject* saveinstrument(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    PyObject* list;
    Int epochType = 0;
    if(!PyArg_ParseTuple(args, "sOi", &s, &list, &epochType))
      throw(Exception("Unable to parse arguments."));
    std::string fname(s);

    UInt arcCount = PyList_Size(list);
    Epoch::Type type = static_cast<Epoch::Type>(epochType);

    std::vector<Arc> arcList(arcCount);
    for(UInt arcNo = 0; arcNo<arcCount; arcNo++)
    {
      Matrix M = fromPyObject(PyList_GetItem(list, arcNo));
      arcList[arcNo] = Arc(M, type);
    }

    InstrumentFile fileInstrument;
    fileInstrument.write(FileName(fname), arcList);

    Py_RETURN_NONE;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Load a GROOPS arc list
 */
static PyObject* loadarclist(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments."));
    std::string fname(s);

    std::vector<UInt> arcsInterval;
    std::vector<Time> timesInterval;

    readFileArcList(FileName(fname), arcsInterval, timesInterval);

    PyObject* list = PyTuple_New(2);
    PyObject* aI = PyTuple_New(arcsInterval.size());
    PyObject* tI = PyTuple_New(timesInterval.size());

    PyTuple_SetItem(list, 0, aI);
    PyTuple_SetItem(list, 1, tI);

    for(UInt k = 0; k<arcsInterval.size(); k++)
    {
      PyTuple_SetItem(aI, k, Py_BuildValue("i", arcsInterval.at(k)));
      PyTuple_SetItem(tI, k, Py_BuildValue("d", timesInterval.at(k).mjd()));
    }

    return list;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Load GROOPS normal equation info
 */
static PyObject* loadnormalsinfo(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments."));

    Matrix n;
    NormalEquationInfo info;
    readFileNormalEquation(FileName(std::string(s)), info, n);

    PyObject* parameter_names = PyTuple_New(info.parameterName.size());
    for(UInt k = 0; k<info.parameterName.size(); k++)
      PyTuple_SetItem(parameter_names, k, Py_BuildValue("s", info.parameterName.at(k).str().c_str()));

    Vector blockIndex(info.blockIndex.size());
    for(UInt k = 0; k<info.blockIndex.size(); k++)
      blockIndex(k) = info.blockIndex.at(k);

    PyObject* list = PyTuple_New(5);
    PyTuple_SetItem(list, 0, fromMatrix(info.lPl));
    PyTuple_SetItem(list, 1, Py_BuildValue("i", info.observationCount));
    PyTuple_SetItem(list, 2, parameter_names);
    PyTuple_SetItem(list, 3, fromMatrix(blockIndex));
    PyTuple_SetItem(list, 4, fromMatrix(info.usedBlocks));

    return list;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Load GROOPS normal equation
 */
static PyObject* loadnormals(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments."));

    Matrix N, n;
    NormalEquationInfo info;
    readFileNormalEquation(FileName(std::string(s)), info, N, n);

    PyObject* list = PyTuple_New(4);
    PyTuple_SetItem(list, 0, fromMatrix(N));
    PyTuple_SetItem(list, 1, fromMatrix(n));
    PyTuple_SetItem(list, 2, fromMatrix(info.lPl));
    PyTuple_SetItem(list, 3, Py_BuildValue("i", info.observationCount));

    return list;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Save matrix and metadata to GROOPS normals file format
 */
static PyObject* savenormals(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    PyObject *matrix;
    PyObject *rhs;
    PyObject *lPl;
    int obsCount = 0;
    if(!PyArg_ParseTuple(args, "sOOOi", &s, &matrix, &rhs, &lPl, &obsCount))
      throw(Exception("Unable to parse arguments."));

    NormalEquationInfo info;

    Matrix N = fromPyObject(matrix, Matrix::SYMMETRIC, Matrix::UPPER);
    Matrix n = fromPyObject(rhs);
    info.lPl = fromPyObject(lPl);
    info.observationCount = obsCount;
    info.parameterName.resize(n.rows());
    info.blockIndex = std::vector<UInt>({0, n.rows()});
    info.usedBlocks = Matrix(1, 1, 1.0);

    writeFileNormalEquation(FileName(std::string(s)), info, N, n);

    Py_RETURN_NONE;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Read a GROOPS polygon list from file.
 */
static PyObject* loadpolygon(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments."));

    std::vector<Polygon> poly;
    readFilePolygon(FileName(std::string(s)), poly);

    PyObject* polygons = PyTuple_New(poly.size());
    for(UInt k = 0; k < poly.size(); k++)
    {
      Matrix data(poly.at(k).L.size(), 2);
      copy(poly.at(k).L, data.column(0));
      copy(poly.at(k).B, data.column(1));
      PyTuple_SetItem(polygons, k, fromMatrix(data));
    }

    return polygons;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Save a GROOPS polygon list to file.
 */
static PyObject* savepolygon(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    PyObject *list;
    if(!PyArg_ParseTuple(args, "sO", &s, &list))
      throw(Exception("Unable to parse arguments."));

    UInt count = PyTuple_Size(list);
    std::vector<Polygon> poly(count);
    for(UInt k = 0; k < poly.size(); k++)
    {
      Matrix M = fromPyObject(PyTuple_GetItem(list, k));
      poly.at(k).L = M.column(0);
      poly.at(k).B = M.column(1);
    }

    writeFilePolygon(FileName(std::string(s)), poly);

    Py_RETURN_NONE;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Read GROOPS parameter names from file.
 */
static PyObject* loadparameternames(PyObject* /*self*/, PyObject* args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments."));

    std::vector<ParameterName> parameterNames;
    readFileParameterName(FileName(std::string(s)), parameterNames);

    PyObject* parameterNameTuple = PyTuple_New(parameterNames.size());
    for(UInt k = 0; k < parameterNames.size(); k++)
    {
      PyObject *str;
      str = PyString_FromStringAndSize(parameterNames.at(k).str().c_str(), parameterNames.at(k).str().size());
      PyTuple_SetItem(parameterNameTuple, k, str);
    }

    return parameterNameTuple;
  }
  catch(std::exception& e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

/*
 * Read groops signalbias file
 */
static PyObject* loadgnsssignalbias(PyObject*, PyObject *args)
{
  try
  {
    const char *s;
    if(!PyArg_ParseTuple(args, "s", &s))
      throw(Exception("Unable to parse arguments"));
    std::string fname(s);
    GnssSignalBias biases;
    readFileGnssSignalBias(fname, biases);
    PyObject* return_tuple = PyTuple_New(biases.biases.size());
    for(UInt i = 0; i < biases.biases.size(); i++)
    {
      PyObject* current_tuple = PyTuple_New(2);
      PyTuple_SetItem(current_tuple, 0, Py_BuildValue("s", biases.types.at(i).str().c_str()));
      PyTuple_SetItem(current_tuple, 1, Py_BuildValue("d", biases.biases.at(i)));

      PyTuple_SetItem(return_tuple, i, current_tuple);
    }
    return return_tuple;
  }
  catch(std::exception &e)
  {
    PyErr_SetString(groopsError, e.what());
    return NULL;
  }
}

#endif
