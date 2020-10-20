/*
 * Python module for GROOPS file I/O.
 *
 * Copyright (C) 2020 Andreas Kvas
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

#include "groopsio.h"

extern "C" {

  /* Exposed functions */
  static PyMethodDef iobase_methods[] = {
    {"loadmat", (PyCFunction)loadmat, METH_VARARGS, ""},
    {"savemat", (PyCFunction)savemat, METH_VARARGS, ""},

    {"loadgridrectangular", (PyCFunction)loadgridrectangular, METH_VARARGS, ""},
    {"loadgrid", (PyCFunction)loadgrid, METH_VARARGS, ""},
    {"savegrid", (PyCFunction)savegrid, METH_VARARGS, ""},

    {"loadgravityfield", (PyCFunction)loadgravityfield, METH_VARARGS, ""},
    {"savegravityfield", (PyCFunction)savegravityfield, METH_VARARGS, ""},

    {"loadtimesplines", (PyCFunction)loadtimesplines, METH_VARARGS, ""},

    {"loadarclist", (PyCFunction)loadarclist, METH_VARARGS, ""},

    {"loadinstrument", (PyCFunction)loadinstrument, METH_VARARGS, ""},
    {"saveinstrument", (PyCFunction)saveinstrument, METH_VARARGS, ""},

    {"loadstarcamera", (PyCFunction)loadstarcamera, METH_VARARGS, ""},

    {"loadnormalsinfo", (PyCFunction)loadnormalsinfo, METH_VARARGS, ""},
    {"loadnormals", (PyCFunction)loadnormals, METH_VARARGS, ""},
    {"loadparameternames", (PyCFunction)loadparameternames, METH_VARARGS, ""},
    {"savenormals", (PyCFunction)savenormals, METH_VARARGS, ""},

    {"loadpolygon", (PyCFunction)loadpolygon, METH_VARARGS, ""},

    {NULL, NULL, 0, NULL}
  };

  /* Module definition */
  static struct PyModuleDef iobase_definition = {
    PyModuleDef_HEAD_INIT, /*m_base*/
    "groopsiobase",        /*m_name*/
    NULL,                  /*m_doc*/
    -1,                    /*m_size*/
    iobase_methods,        /*m_methods*/
    NULL,                  /*m_slots*/
    NULL,                  /*m_traverse*/
    NULL,                  /*m_clear*/
    NULL                   /*m_free*/
  };

  /* Initialization function for the module */
  PyMODINIT_FUNC PyInit_groopsiobase()
  {
    PyObject *module = PyModule_Create(&iobase_definition);
    import_array();

    groopsError = PyErr_NewException("groopsiobase.GROOPSError", NULL, NULL);
    Py_INCREF(groopsError);
    PyModule_AddObject(module, "GROOPSError", groopsError);

    return module;
  }
}
