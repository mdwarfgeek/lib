#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>

#include <stdlib.h>
#include <inttypes.h>
#include <float.h>

#include "ap.h"
#include "cvtunit.h"
#include "lfa.h"

/* Wrappers for the observer and source structures */

struct lfa_source_object {
  PyObject_HEAD
  struct source s;
};

struct lfa_observer_object {
  PyObject_HEAD
  struct observer o;
  struct dtai_table dtab;
  struct iers_table itab;
  struct jpleph_table jtab;
  struct jpleph_table ttab;
  struct jpleph_table *tptr;
};

struct lfa_ap_object {
  PyObject_HEAD
  int allocated;
  struct ap ap;
};

static int lfa_source_init (struct lfa_source_object *self,
                            PyObject *args,
                            PyObject *kwds) {
  static char *kwlist[] = { "ra", "de",
                            "pmra", "pmde",
                            "plx", "vrad",
                            "epoch",
                            NULL };

  int narg, iplan;

  double ra, de;
  double pmra = 0, pmde = 0;
  double plx = 0, vrad = 0;
  double epoch = 2000.0;

  /* How many non-keyword arguments? */
  narg = PyTuple_Size(args);

  if(narg == 1) {
    if(!PyArg_ParseTuple(args, "i", &iplan))
      return(-1);

    self->s.type = iplan;
  }
  else {
    if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                    "dd|ddddd", kwlist,
                                    &ra, &de,
                                    &pmra, &pmde,
                                    &plx, &vrad,
                                    &epoch))
      return(-1);
    
    source_star(&(self->s),
                ra, de,
                pmra, pmde,
                plx, vrad,
                epoch);
  }

  return(0);
}

static PyObject *lfa_source_ecliptic (struct lfa_source_object *self,
                                      PyObject *args,
                                      PyObject *kwds) {
  double sec[3], l, b;

  /* NOTE: coordinates are mean ecliptic of J2000 at catalogue epoch,
     which is what I'm assuming most people mean when they say "ecliptic
     coordinates".  BUT NOT ALWAYS! */
  m_x_v(gcrs2ecl, self->s.ref_n, sec);

  v_to_ad(sec, 0, &l, &b);
  l = ranormp(l);

  return(Py_BuildValue("dd", l, b));
}

static PyObject *lfa_source_galactic (struct lfa_source_object *self,
                                      PyObject *args,
                                      PyObject *kwds) {
  PyObject *uvwarr = NULL;
  npy_intp outdim[1] = { 3 };

  double sgal[3], l, b, *uvw;
  int i;

  uvwarr = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!uvwarr)
    goto error;

  uvw = (double *) PyArray_DATA(uvwarr);

  /* At catalogue epoch */
  m_x_v(eq2gal, self->s.ref_n, sgal);
  m_x_v(eq2gal, self->s.ref_dndt, uvw);

  v_to_ad(sgal, 0, &l, &b);
  l = ranormp(l);

  for(i = 0; i < 3; i++)
    uvw[i] *= AU / (DAY*1000*self->s.pr);  /* to km/s */

  return(Py_BuildValue("ddN", l, b,
                       PyArray_Return((PyArrayObject *) uvwarr)));

 error:
  PyArray_XDECREF_ERR((PyArrayObject *) uvwarr);

  return(NULL);
}

static PyMemberDef lfa_source_members[] = {
  { NULL, 0, 0, 0, NULL }
};

static PyMethodDef lfa_source_methods[] = {
  { "galactic", (PyCFunction) lfa_source_galactic,
    METH_VARARGS | METH_KEYWORDS,
    "(l, b, uvw) = galactic()" },
  { "ecliptic", (PyCFunction) lfa_source_ecliptic,
    METH_VARARGS | METH_KEYWORDS,
    "(l, b) = ecliptic()" },
  { NULL, NULL, 0, NULL }
};

static PyTypeObject lfa_source_type = {
  PyObject_HEAD_INIT(NULL)
  0,                                   /* ob_size */
  "lfa.source",                        /* tp_name */
  sizeof(struct lfa_source_object),    /* tp_basicsize */
  0,                                   /* tp_itemsize */
  0,                                   /* tp_dealloc */
  0,                                   /* tp_print */
  0,                                   /* tp_getattr */
  0,                                   /* tp_setattr */
  0,                                   /* tp_compare */
  0,                                   /* tp_repr */
  0,                                   /* tp_as_number */
  0,                                   /* tp_as_sequence */
  0,                                   /* tp_as_mapping */
  0,                                   /* tp_hash */
  0,                                   /* tp_call */
  0,                                   /* tp_str */
  0,                                   /* tp_getattro */
  0,                                   /* tp_setattro */
  0,                                   /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,  /* tp_flags */
  "Source structure",                  /* tp_doc */
  0,                                   /* tp_traverse */
  0,                                   /* tp_clear */
  0,                                   /* tp_richcompare */
  0,                                   /* tp_weaklistoffset */
  0,                                   /* tp_iter */
  0,                                   /* tp_iternext */
  lfa_source_methods,                  /* tp_methods */
  lfa_source_members,                  /* tp_members */
  0,                                   /* tp_getset */
  0,                                   /* tp_base */
  0,                                   /* tp_dict */
  0,                                   /* tp_descr_get */
  0,                                   /* tp_descr_set */
  0,                                   /* tp_dictoffset */
  (initproc) lfa_source_init,          /* tp_init */
  0,                                   /* tp_alloc */
  NULL                                 /* tp_new */
};

static void lfa_observer_error (int rv, char *func, char *env, char *file) {
  char buf[1024];

  switch(rv) {
  case -1:
    if(file)
      snprintf(buf, sizeof(buf),
               "%s: %s: %s", func, file, strerror(errno));
    else
      snprintf(buf, sizeof(buf),
               "%s: %s", func, strerror(errno));
    break;
  case -2:
    snprintf(buf, sizeof(buf),
             "%s: %s environment variable not set", func, env);
    break;
  case -3:
    if(file)
      snprintf(buf, sizeof(buf),
               "%s: parse error, %s may be corrupted", func, file);
    else
      snprintf(buf, sizeof(buf),
               "%s: parse error, file may be corrupted", func);
    break;
  default:
    snprintf(buf, sizeof(buf),
             "%s failed", func);
    break;
  }

  PyErr_SetString(PyExc_RuntimeError, buf);
}

static PyObject *lfa_observer_new (PyTypeObject *type,
                                   PyObject *args,
                                   PyObject *kwds) {
  struct lfa_observer_object *self = NULL;
  int rv;

  /* Allocate */
  self = (struct lfa_observer_object *) type->tp_alloc(type, 0);
  if(self) {
    /* Setup Earth orientation data and JPL ephemerides */
    rv = dtai_read(&(self->dtab), (char *) NULL);
    if(rv) {
      lfa_observer_error(rv, "dtai_read", "IERS_DATA", "tai-utc.dat");
      goto error;
    }

    rv = iers_open(&(self->itab), &(self->dtab), (char *) NULL);
    if(rv) {
      lfa_observer_error(rv, "iers_open", "IERS_DATA", "finals2000A.data");
      goto error;
    }

    rv = jpleph_open(&(self->jtab), 0, (char *) NULL);
    if(rv) {
      lfa_observer_error(rv, "jpleph_open", "JPLEPH_DATA", (char *) NULL);
      goto error;
    }

    if(!self->jtab.has_time) {
      rv = jpleph_open(&(self->ttab), 1, (char *) NULL);
      if(rv) {
        lfa_observer_error(rv, "jpleph_open", "TIMEEPH_DATA", (char *) NULL);
        goto error;
      }
      
      self->tptr = &(self->ttab);
    }
    else
      self->tptr = NULL;
  }

  return((PyObject *) self);

 error:

  return(NULL);
}

static int lfa_observer_init (struct lfa_observer_object *self,
                              PyObject *args,
                              PyObject *kwds) {
  static char *kwlist[] = { "longitude",
                            "latitude",
                            "height",
                            NULL };
  double longitude, latitude, height = 0;

  if(PyTuple_Size(args)) {
    /* Get arguments */
    if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                    "dd|d", kwlist,
                                    &longitude,
                                    &latitude,
                                    &height))
      return(-1);
    
    observer_init(&(self->o), longitude, latitude, height);
  }
  else
    observer_geoc(&(self->o));

  return(0);
}

static void lfa_observer_dealloc (struct lfa_observer_object *self) {
  
  if(self->tptr)
    jpleph_close(self->tptr);

  jpleph_close(&(self->jtab));
  iers_close(&(self->itab));
  dtai_free(&(self->dtab));

  self->ob_type->tp_free((PyObject *) self);
}

static PyObject *lfa_observer_ast2obs (struct lfa_observer_object *self,
                                       PyObject *args,
                                       PyObject *kwds) {
  static char *kwlist[] = { "s",
                            "dsdt",
                            "pr",
                            "mask",
                            NULL };
  PyObject *sarg, *dsdtarg = NULL;
  double pr = 0;
  int mask = TR_TO_OBS_AZ;

  PyObject *sarr = NULL, *dsdtarr = NULL;
  PyObject *sout = NULL, *dsdtout = NULL;

  npy_intp outdim[1] = { 3 };

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O|Odi", kwlist,
                                  &sarg,
                                  &dsdtarg,
                                  &pr,
                                  &mask))
    goto error;

  sarr = PyArray_FROM_OTF(sarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!sarr)
    goto error;

  if(PyArray_Size(sarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 's' is too small");
    goto error;
  }

  sout = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!sout)
    goto error;

  memcpy(PyArray_DATA(sout), PyArray_DATA(sarr), 3*sizeof(double));

  if(dsdtarg && dsdtarg != Py_None) {
    dsdtarr = PyArray_FROM_OTF(dsdtarg, NPY_DOUBLE,
                               NPY_IN_ARRAY | NPY_FORCECAST);
    if(!dsdtarr)
      goto error;

    if(PyArray_Size(dsdtarr) < 3) {
      PyErr_SetString(PyExc_IndexError,
                      "array 'dsdt' is too small");
      goto error;
    }

    dsdtout = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
    if(!dsdtout)
      goto error;
    
    memcpy(PyArray_DATA(dsdtout), PyArray_DATA(dsdtarr), 3*sizeof(double));
  }

  observer_ast2obs(&(self->o),
                   PyArray_DATA(sout),
                   dsdtout ? PyArray_DATA(dsdtout) : NULL,
                   pr,
                   mask);

  Py_DECREF(sarr);
  Py_XDECREF(dsdtarr);

  if(dsdtout)
    return(Py_BuildValue("NN",
                         PyArray_Return((PyArrayObject *) sout),
                         PyArray_Return((PyArrayObject *) dsdtout)));
  else
    return(Py_BuildValue("N",
                         PyArray_Return((PyArrayObject *) sout)));

 error:
  Py_XDECREF(sarr);
  Py_XDECREF(dsdtarr);

  PyArray_XDECREF_ERR((PyArrayObject *) sout);
  PyArray_XDECREF_ERR((PyArrayObject *) dsdtout);

  return(NULL);
}

static PyObject *lfa_observer_bary_delay (struct lfa_observer_object *self,
                                          PyObject *args,
                                          PyObject *kwds) {
  static char *kwlist[] = { "s",
                            "pr",
                            NULL };
  PyObject *sarg;
  double pr = 0;

  PyObject *sarr = NULL;
  double rv;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O|d", kwlist,
                                  &sarg,
                                  &pr))
    goto error;

  sarr = PyArray_FROM_OTF(sarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!sarr)
    goto error;

  if(PyArray_Size(sarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 's' is too small");
    goto error;
  }

  rv = bary_delay(&(self->o), PyArray_DATA(sarr), pr);

  Py_DECREF(sarr);

  return(Py_BuildValue("d", rv));

 error:
  Py_XDECREF(sarr);

  return(NULL);
}

static PyObject *lfa_observer_bary_doppler (struct lfa_observer_object *self,
                                            PyObject *args,
                                            PyObject *kwds) {
  static char *kwlist[] = { "src",
                            "s",
                            "dsdt",
                            "pr",
                            NULL };
  PyObject *src_arg, *sarg, *dsdtarg;
  double pr = 0;

  struct lfa_source_object *src;
  PyObject *sarr = NULL, *dsdtarr = NULL;
  double rv;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OOO|d", kwlist,
                                  &src_arg, &sarg, &dsdtarg,
                                  &pr))
    goto error;

  if(!PyObject_TypeCheck(src_arg, &lfa_source_type)) {
    PyErr_SetString(PyExc_TypeError, "argument 1 must be an lfa.source");
    goto error;
  }

  src = (struct lfa_source_object *) src_arg;

  sarr = PyArray_FROM_OTF(sarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!sarr)
    goto error;

  if(PyArray_Size(sarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 's' is too small");
    goto error;
  }

  dsdtarr = PyArray_FROM_OTF(dsdtarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!dsdtarr)
    goto error;

  if(PyArray_Size(dsdtarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'dsdt' is too small");
    goto error;
  }

  rv = bary_doppler(&(self->o), src->s.ref_n,
                    PyArray_DATA(sarr), PyArray_DATA(dsdtarr), pr);

  Py_DECREF(sarr);
  Py_DECREF(dsdtarr);

  return(Py_BuildValue("d", rv));

 error:
  Py_XDECREF(sarr);
  Py_XDECREF(dsdtarr);

  return(NULL);
}

static PyObject *lfa_observer_dtai (struct lfa_observer_object *self,
                                    PyObject *args,
                                    PyObject *kwds) {
  static char *kwlist[] = { "mjd",
                            "frac",
                            NULL };
  int mjd;
  double frac, rv;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "id", kwlist,
                                  &mjd,
                                  &frac))
    return(NULL);

  rv = dtai(self->itab.dtai_tab, mjd, frac);

  return(Py_BuildValue("d", rv));
}

static PyObject *lfa_observer_obs2ast (struct lfa_observer_object *self,
                                       PyObject *args,
                                       PyObject *kwds) {
  static char *kwlist[] = { "s",
                            "pr",
                            "mask",
                            NULL };
  PyObject *sarg;
  double pr = 0;
  int mask = TR_TO_OBS_AZ;

  PyObject *sarr = NULL;
  PyObject *sout = NULL;

  npy_intp outdim[1] = { 3 };

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O|di", kwlist,
                                  &sarg,
                                  &pr,
                                  &mask))
    goto error;

  sarr = PyArray_FROM_OTF(sarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!sarr)
    goto error;

  if(PyArray_Size(sarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 's' is too small");
    goto error;
  }

  sout = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!sout)
    goto error;

  memcpy(PyArray_DATA(sout), PyArray_DATA(sarr), 3*sizeof(double));

  observer_obs2ast(&(self->o),
                   PyArray_DATA(sout),
                   pr,
                   mask);

  Py_DECREF(sarr);

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) sout)));

 error:
  Py_XDECREF(sarr);
  PyArray_XDECREF_ERR((PyArrayObject *) sout);

  return(NULL);
}

static PyObject *lfa_observer_place (struct lfa_observer_object *self,
                                     PyObject *args,
                                     PyObject *kwds) {
  static char *kwlist[] = { "src",
                            "mask",
                            NULL };
  PyObject *src_arg;
  struct lfa_source_object *src;
  int mask;

  npy_intp outdim[1] = { 3 };
  PyObject *s_ret = NULL, *dsdt_ret = NULL;
  double pr;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "Oi", kwlist,
                                  &src_arg,
                                  &mask))
    goto error;

  if(!PyObject_TypeCheck(src_arg, &lfa_source_type)) {
    PyErr_SetString(PyExc_TypeError, "argument 1 must be an lfa.source");
    goto error;
  }

  src = (struct lfa_source_object *) src_arg;

  s_ret = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!s_ret)
    goto error;

  dsdt_ret = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!dsdt_ret)
    goto error;

  source_place(&(self->o), &(src->s), &(self->jtab),
               (unsigned char) mask,
               PyArray_DATA(s_ret),
               PyArray_DATA(dsdt_ret),
               &pr);

  return(Py_BuildValue("NNd",
                       PyArray_Return((PyArrayObject *) s_ret),
                       PyArray_Return((PyArrayObject *) dsdt_ret),
                       pr));

 error:
  PyArray_XDECREF_ERR((PyArrayObject *) s_ret);
  PyArray_XDECREF_ERR((PyArrayObject *) dsdt_ret);

  return(NULL);
}

static PyObject *lfa_observer_refract (struct lfa_observer_object *self,
                                       PyObject *args,
                                       PyObject *kwds) {
  static char *kwlist[] = { "temperat",
                            "humidity",
                            "pressure",
                            "wavelength",
                            NULL };
  double temperat, humidity, pressure, wavelength;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "dddd", kwlist,
                                  &temperat,
                                  &humidity,
                                  &pressure,
                                  &wavelength))
    return(NULL);

  refract_const(temperat, humidity, pressure, wavelength,
                self->o.height, self->o.refco);


  return(Py_None);
}

static PyObject *lfa_observer_update (struct lfa_observer_object *self,
                                      PyObject *args,
                                      PyObject *kwds) {
  static char *kwlist[] = { "utc",
                            "ttmutc",
                            "mask",
                            NULL };
  double utc, ttmutc;
  int mask = OBSERVER_UPDATE_ALL;
  int rv;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "dd|i", kwlist,
                                  &utc,
                                  &ttmutc,
                                  &mask))
    return(NULL);

  rv = observer_update(&(self->o),
                       &(self->jtab),
                       self->tptr,
                       &(self->itab),
                       utc, ttmutc, (unsigned char) mask);

  return(Py_BuildValue("i", rv));
}

static PyMemberDef lfa_observer_members[] = {
  { "latitude",  T_DOUBLE, offsetof(struct lfa_observer_object, o.latitude),   0, "latitude" },
  { "longitude", T_DOUBLE, offsetof(struct lfa_observer_object, o.longitude),  0, "longitude" },
  { "height",    T_DOUBLE, offsetof(struct lfa_observer_object, o.height),     0, "height" },
  { "utc",  T_DOUBLE, offsetof(struct lfa_observer_object, o.utc),  0, "utc" },
  { "tt",   T_DOUBLE, offsetof(struct lfa_observer_object, o.tt),   0, "tt" },
  { "tdb",  T_DOUBLE, offsetof(struct lfa_observer_object, o.tdb),  0, "tdb" },
  { "dtdb", T_DOUBLE, offsetof(struct lfa_observer_object, o.dtdb), 0, "dtdb" },
  { NULL, 0, 0, 0, NULL }
};

static PyMethodDef lfa_observer_methods[] = {
  { "ast2obs", (PyCFunction) lfa_observer_ast2obs,
    METH_VARARGS | METH_KEYWORDS,
    "(s, dsdt) = ast2obs(s, dsdt=None, pr=0, mask=lfa.TR_TO_OBS_AZ)" },
  { "bary_delay", (PyCFunction) lfa_observer_bary_delay,
    METH_VARARGS | METH_KEYWORDS,
    "dly = bary_delay(s, pr=0)" },
  { "bary_doppler", (PyCFunction) lfa_observer_bary_doppler,
    METH_VARARGS | METH_KEYWORDS,
    "shift = bary_doppler(src, s=None, dsdt=None, pr=0)" },
  { "dtai", (PyCFunction) lfa_observer_dtai,
    METH_VARARGS | METH_KEYWORDS,
    "ttmutc = dtai(mjd, frac)" },
  { "obs2ast", (PyCFunction) lfa_observer_obs2ast,
    METH_VARARGS | METH_KEYWORDS,
    "s = obs2ast(s, pr=0, mask=lfa.TR_TO_OBS_AZ)" },
  { "place", (PyCFunction) lfa_observer_place,
    METH_VARARGS | METH_KEYWORDS,
    "(s, dsdt, pr) = place(src, mask)" },
  { "refract", (PyCFunction) lfa_observer_refract,
    METH_VARARGS | METH_KEYWORDS,
    "refract(temperat, humidity, pressure, wavelength)" },
  { "update", (PyCFunction) lfa_observer_update,
    METH_VARARGS | METH_KEYWORDS,
    "update(utc, ttmutc, mask=lfa.OBSERVER_UPDATE_ALL)" },
  { NULL, NULL, 0, NULL }
};

static PyTypeObject lfa_observer_type = {
  PyObject_HEAD_INIT(NULL)
  0,                                   /* ob_size */
  "lfa.observer",                      /* tp_name */
  sizeof(struct lfa_observer_object),  /* tp_basicsize */
  0,                                   /* tp_itemsize */
  (destructor) lfa_observer_dealloc,   /* tp_dealloc */
  0,                                   /* tp_print */
  0,                                   /* tp_getattr */
  0,                                   /* tp_setattr */
  0,                                   /* tp_compare */
  0,                                   /* tp_repr */
  0,                                   /* tp_as_number */
  0,                                   /* tp_as_sequence */
  0,                                   /* tp_as_mapping */
  0,                                   /* tp_hash */
  0,                                   /* tp_call */
  0,                                   /* tp_str */
  0,                                   /* tp_getattro */
  0,                                   /* tp_setattro */
  0,                                   /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,  /* tp_flags */
  "Observer structure",                /* tp_doc */
  0,                                   /* tp_traverse */
  0,                                   /* tp_clear */
  0,                                   /* tp_richcompare */
  0,                                   /* tp_weaklistoffset */
  0,                                   /* tp_iter */
  0,                                   /* tp_iternext */
  lfa_observer_methods,                /* tp_methods */
  lfa_observer_members,                /* tp_members */
  0,                                   /* tp_getset */
  0,                                   /* tp_base */
  0,                                   /* tp_dict */
  0,                                   /* tp_descr_get */
  0,                                   /* tp_descr_set */
  0,                                   /* tp_dictoffset */
  (initproc) lfa_observer_init,        /* tp_init */
  0,                                   /* tp_alloc */
  lfa_observer_new                     /* tp_new */
};

static PyObject *lfa_ap_new (PyTypeObject *type,
                             PyObject *args,
                             PyObject *kwds) {
  struct lfa_ap_object *self = NULL;

  /* Allocate */
  self = (struct lfa_ap_object *) type->tp_alloc(type, 0);
  if(self)
    self->allocated = 0;

  return((PyObject *) self);
}

static void lfa_ap_dealloc (struct lfa_ap_object *self) {
  if(self->allocated) {
    ap_free(&(self->ap));
    self->allocated = 0;
  }

  self->ob_type->tp_free((PyObject *) self);
}

static int lfa_ap_init (struct lfa_ap_object *self,
                        PyObject *args,
                        PyObject *kwds) {
  static char *kwlist[] = { "nx",
                            "ny",
                            NULL };
  int nx, ny;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "ii", kwlist,
                                  &nx, &ny))
    return(-1);
    
  if(ap_init(&(self->ap), nx, ny)) {
    PyErr_SetString(PyExc_MemoryError, "ap_init");
    return(-1);
  }

  self->allocated = 1;

  return(0);
}

static PyObject *lfa_ap_image (struct lfa_ap_object *self,
                               PyObject *args,
                               PyObject *kwds) {
  static char *kwlist[] = { "map",
                            "minpix",
                            "sky",
                            "thresh",
                            "filtmap",
                            "mask",
                            NULL };
  PyObject *maparg;
  PyObject *filtmaparg = NULL;
  PyObject *maskarg = NULL;
  int minpix;
  float sky, thresh;

  PyObject *maparr = NULL;
  PyObject *filtmaparr = NULL;
  PyObject *maskarr = NULL;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "Oiff|OO", kwlist,
                                  &maparg,
                                  &minpix,
                                  &sky,
                                  &thresh,
                                  &filtmaparg,
                                  &maskarg))
    return(NULL);

  maparr = PyArray_FROM_OTF(maparg, NPY_FLOAT, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!maparr)
    goto error;

  if(PyArray_Size(maparr) < self->ap.nx * self->ap.ny) {
    PyErr_SetString(PyExc_IndexError, "map too small");
    goto error;
  }

  if(filtmaparg) {
    filtmaparr = PyArray_FROM_OTF(filtmaparg, NPY_FLOAT, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!filtmaparr)
      goto error;

    if(PyArray_Size(filtmaparr) < self->ap.nx * self->ap.ny) {
      PyErr_SetString(PyExc_IndexError, "filtmap too small");
      goto error;
    }
  }

  if(maskarg) {
    maskarr = PyArray_FROM_OTF(maskarg, NPY_UINT8, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!maskarr)
      goto error;

    if(PyArray_Size(maskarr) < self->ap.nx * self->ap.ny) {
      PyErr_SetString(PyExc_IndexError, "mask too small");
      goto error;
    }
  }

  if(ap_image(&(self->ap),
              PyArray_DATA(maparr),
              filtmaparr ? PyArray_DATA(filtmaparr) : NULL,
              maskarr ? PyArray_DATA(maskarr) : NULL,
              minpix, sky, thresh)) {
    PyErr_SetString(PyExc_RuntimeError, "ap_image");
    goto error;
  }

  Py_DECREF(maparr);
  Py_XDECREF(filtmaparr);
  Py_XDECREF(maskarr);

  return(Py_BuildValue("i", self->ap.nsources));

 error:
  Py_XDECREF(maparr);
  Py_XDECREF(filtmaparr);
  Py_XDECREF(maskarr);

  return(NULL);
}

static PyObject *lfa_ap_source (struct lfa_ap_object *self,
                                PyObject *args,
                                PyObject *kwds) {
  static char *kwlist[] = { "isrc",
                            NULL };

  int isrc;
  struct ap_source *src;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "i", kwlist,
                                  &isrc))
    return(NULL);

  if(!self->allocated ||
     isrc < 0 ||
     isrc >= self->ap.nsources) {
    PyErr_SetString(PyExc_IndexError, "isrc out of range");
    goto error;
  }

  src = self->ap.sources + isrc;

  return(Py_BuildValue("ddddi", src->x, src->y, src->ziso, src->zpk, src->areal[0]));

 error:

  return(NULL);
}

static PyObject *lfa_ap_ellipse (struct lfa_ap_object *self,
                                 PyObject *args,
                                 PyObject *kwds) {
  static char *kwlist[] = { "isrc",
                            NULL };

  int isrc;
  struct ap_source *src;
  double gfwhm, ell, spa, cpa;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "i", kwlist,
                                  &isrc))
    return(NULL);

  if(!self->allocated ||
     isrc < 0 ||
     isrc >= self->ap.nsources) {
    PyErr_SetString(PyExc_IndexError, "isrc out of range");
    goto error;
  }

  src = self->ap.sources + isrc;

  ap_ellipse(src, &gfwhm, &ell, &spa, &cpa);

  return(Py_BuildValue("dddd", gfwhm, ell, spa, cpa));

 error:

  return(NULL);
}

static PyMemberDef lfa_ap_members[] = {
  { NULL, 0, 0, 0, NULL }
};

static PyMethodDef lfa_ap_methods[] = {
  { "image", (PyCFunction) lfa_ap_image,
    METH_VARARGS | METH_KEYWORDS,
    "nsources = image(map, minpix, sky, thresh)" },
  { "source", (PyCFunction) lfa_ap_source,
    METH_VARARGS | METH_KEYWORDS,
    "x, y, ziso, zpk, niso = source(isrc)" },
  { "ellipse", (PyCFunction) lfa_ap_ellipse,
    METH_VARARGS | METH_KEYWORDS,
    "gfwhm, ell, spa, cps = ellipse(isrc)" },
  { NULL, NULL, 0, NULL }
};

static PyTypeObject lfa_ap_type = {
  PyObject_HEAD_INIT(NULL)
  0,                                   /* ob_size */
  "lfa.ap",                            /* tp_name */
  sizeof(struct lfa_ap_object),        /* tp_basicsize */
  0,                                   /* tp_itemsize */
  (destructor) lfa_ap_dealloc,         /* tp_dealloc */
  0,                                   /* tp_print */
  0,                                   /* tp_getattr */
  0,                                   /* tp_setattr */
  0,                                   /* tp_compare */
  0,                                   /* tp_repr */
  0,                                   /* tp_as_number */
  0,                                   /* tp_as_sequence */
  0,                                   /* tp_as_mapping */
  0,                                   /* tp_hash */
  0,                                   /* tp_call */
  0,                                   /* tp_str */
  0,                                   /* tp_getattro */
  0,                                   /* tp_setattro */
  0,                                   /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,  /* tp_flags */
  "ap structure",                      /* tp_doc */
  0,                                   /* tp_traverse */
  0,                                   /* tp_clear */
  0,                                   /* tp_richcompare */
  0,                                   /* tp_weaklistoffset */
  0,                                   /* tp_iter */
  0,                                   /* tp_iternext */
  lfa_ap_methods,                      /* tp_methods */
  lfa_ap_members,                      /* tp_members */
  0,                                   /* tp_getset */
  0,                                   /* tp_base */
  0,                                   /* tp_dict */
  0,                                   /* tp_descr_get */
  0,                                   /* tp_descr_set */
  0,                                   /* tp_dictoffset */
  (initproc) lfa_ap_init,              /* tp_init */
  0,                                   /* tp_alloc */
  lfa_ap_new                           /* tp_new */
};

static PyObject *lfa_source_star_vec (PyObject *self,
                                      PyObject *args,
                                      PyObject *kwds) {
  static char *kwlist[] = { "n",
                            "dndt",
                            "pr",
                            "epoch",
                            NULL };

  PyObject *narg, *dndtarg;
  double pr = 0;
  double epoch = 2000.0;

  PyObject *narr = NULL, *dndtarr = NULL;

  struct lfa_source_object *src = NULL;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OO|dd", kwlist,
                                  &narg, &dndtarg,
                                  &pr, &epoch))
    goto error;

  narr = PyArray_FROM_OTF(narg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!narr)
    goto error;

  if(PyArray_Size(narr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'n' is too small");
    goto error;
  }

  dndtarr = PyArray_FROM_OTF(dndtarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!dndtarr)
    goto error;

  if(PyArray_Size(dndtarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'dndt' is too small");
    goto error;
  }
     
  /* Allocate */
  src = (struct lfa_source_object *)
    lfa_source_type.tp_alloc(&lfa_source_type, 0);
  if(!src)
    goto error;

  source_star_vec(&(src->s),
                  PyArray_DATA(narr),
                  PyArray_DATA(dndtarr),
                  pr, epoch);

  Py_DECREF(narr);
  Py_DECREF(dndtarr);

  return((PyObject *) src);

 error:
  Py_XDECREF(narr);
  Py_XDECREF(dndtarr);
  Py_XDECREF(src);

  return(NULL);
}

static PyObject *lfa_source_elem (PyObject *self,
                                  PyObject *args,
                                  PyObject *kwds) {
  static char *kwlist[] = { "eltype",
                            "epoch", "incl", "anode", "longperi",
                            "aq", "ecc", "lm", "nn",
                            NULL };

  struct lfa_source_object *src = NULL;

  unsigned char eltype;
  double epoch, incl, anode, longperi, aq, ecc, lm, nn;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "bddddddddd", kwlist,
                                  &eltype,
                                  &epoch, &incl, &anode, &longperi,
                                  &aq, &ecc, &lm, &nn))
    return(NULL);

  /* Allocate */
  src = (struct lfa_source_object *)
    lfa_source_type.tp_alloc(&lfa_source_type, 0);

  if(src) {
    source_elem(&(src->s),
                eltype,
                epoch, incl, anode, longperi,
                aq, ecc, lm, nn);
  }

  return((PyObject *) src);
}

static PyObject *lfa_source_mpc (PyObject *self,
                                 PyObject *args,
                                 PyObject *kwds) {
  static char *kwlist[] = { "line",
                            NULL };
  char *line;

  struct lfa_source_object *src = NULL;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "s", kwlist,
                                  &line))
    return(NULL);

  /* Allocate */
  src = (struct lfa_source_object *)
    lfa_source_type.tp_alloc(&lfa_source_type, 0);
  if(!src)
    goto error;

  if(mpc_convert(line, &(src->s))) {
    PyErr_SetString(PyExc_RuntimeError,
                    "parse error");
  }

  return((PyObject *) src);

 error:
  Py_XDECREF(src);

  return(NULL);
}

static PyObject *lfa_source_planet (PyObject *self,
                                    PyObject *args,
                                    PyObject *kwds) {
  static char *kwlist[] = { "n",
                            NULL };

  struct lfa_source_object *src = NULL;

  unsigned char n;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "b", kwlist,
                                  &n))
    return(NULL);

  /* Allocate */
  src = (struct lfa_source_object *)
    lfa_source_type.tp_alloc(&lfa_source_type, 0);

  if(src) {
    src->s.type = n;
  }

  return((PyObject *) src);
}

static PyObject *lfa_base60_to_10 (PyObject *self,
                                   PyObject *args,
                                   PyObject *kwds) {
  static char *kwlist[] = { "inangle",
                            "sep",
                            "inunit",
                            "outunit",
                            NULL };
  char *inangle, *sep;
  short inunit, outunit;

  int rv;
  double outangle;
  char *ep;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "sshh", kwlist,
                                  &inangle, &sep, &inunit, &outunit))
    goto error;

  if(base60_to_10(inangle, &ep, sep, inunit, &outangle, outunit))
    rv = -1;
  else
    rv = (int) (ep - inangle);

  return(Py_BuildValue("di", outangle, rv));

 error:
  return(NULL);
}

static PyObject *lfa_base10_to_60 (PyObject *self,
                                   PyObject *args,
                                   PyObject *kwds) {
  static char *kwlist[] = { "inangle",
                            "inunit",
                            "sep",
                            "sign",
                            "dp",
                            "outunit",
                            NULL };
  double inangle;
  short inunit, dp, outunit;
  char *sep, *sign;

  char outangle[128];
  int rv;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "dhsshh", kwlist,
                                  &inangle, &inunit,
                                  &sep, &sign, &dp, &outunit))
    goto error;

  rv = base10_to_60(inangle, inunit, outangle, sizeof(outangle),
                    sep, sign, dp, outunit);
  if(rv == -1) {  /* XXX - handle this better */
    PyErr_SetString(PyExc_MemoryError,
                    "ran out of buffer space");
    goto error;
  }

  return(Py_BuildValue("s", outangle));

 error:
  return(NULL);
}

static PyObject *lfa_vec2tp (PyObject *self,
                             PyObject *args,
                             PyObject *kwds) {
  static char *kwlist[] = { "s",
                            "tp",
                            NULL };
  PyObject *sarg, *tparg;

  PyObject *sarr = NULL, *tparr = NULL;
  double x, y;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OO", kwlist,
                                  &sarg, &tparg))
    goto error;

  sarr = PyArray_FROM_OTF(sarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!sarr)
    goto error;

  if(PyArray_Size(sarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 's' is too small");
    goto error;
  }

  tparr = PyArray_FROM_OTF(tparg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!tparr)
    goto error;

  if(PyArray_Size(tparr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'tp' is too small");
    goto error;
  }

  vec2tp(PyArray_DATA(sarr), PyArray_DATA(tparr), &x, &y);

  Py_DECREF(sarr);
  Py_DECREF(tparr);

  return(Py_BuildValue("dd", x, y));

 error:
  Py_XDECREF(sarr);
  Py_XDECREF(tparr);

  return(NULL);
}

static PyObject *lfa_tp2vec (PyObject *self,
                             PyObject *args,
                             PyObject *kwds) {
  static char *kwlist[] = { "x",
                            "y",
                            "tp",
                            NULL };
  PyObject *tparg;

  PyObject *tparr = NULL;
  double x, y;

  PyObject *sout = NULL;

  npy_intp outdim[1] = { 3 };

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "ddO", kwlist,
                                  &x, &y, &tparg))
    goto error;

  tparr = PyArray_FROM_OTF(tparg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!tparr)
    goto error;

  if(PyArray_Size(tparr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'tp' is too small");
    goto error;
  }

  sout = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!sout)
    goto error;

  tp2vec(x, y, PyArray_DATA(tparr), PyArray_DATA(sout));

  Py_DECREF(tparr);

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) sout)));

 error:
  Py_XDECREF(tparr);
  PyArray_XDECREF_ERR((PyArrayObject *) sout);

  return(NULL);
}

static PyObject *lfa_v_airmass (PyObject *self,
                                PyObject *args,
                                PyObject *kwds) {
  static char *kwlist[] = { "s",
                            NULL };
  PyObject *sarg;

  PyObject *sarr = NULL;
  double rv;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O", kwlist,
                                  &sarg))
    goto error;

  sarr = PyArray_FROM_OTF(sarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!sarr)
    goto error;

  if(PyArray_Size(sarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 's' is too small");
    goto error;
  }

  rv = v_airmass(PyArray_DATA(sarr));

  Py_DECREF(sarr);

  return(Py_BuildValue("d", rv));

 error:
  Py_XDECREF(sarr);

  return(NULL);
}

static PyObject *lfa_v_parallactic (PyObject *self,
                                    PyObject *args,
                                    PyObject *kwds) {
  static char *kwlist[] = { "sinphi",
                            "cosphi",
                            "s",
                            NULL };
  double sinphi, cosphi;
  PyObject *sarg;

  PyObject *sarr = NULL;
  double rv;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "ddO", kwlist,
                                  &sinphi, &cosphi, &sarg))
    goto error;

  sarr = PyArray_FROM_OTF(sarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!sarr)
    goto error;

  if(PyArray_Size(sarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 's' is too small");
    goto error;
  }

  rv = v_parallactic(sinphi, cosphi, PyArray_DATA(sarr));

  Py_DECREF(sarr);

  return(Py_BuildValue("d", rv));

 error:
  Py_XDECREF(sarr);

  return(NULL);
}

static PyObject *lfa_ad_to_v (PyObject *self,
                              PyObject *args,
                              PyObject *kwds) {
  static char *kwlist[] = { "a",
                            "d",
                            NULL };
  double a, d;

  PyObject *sout = NULL;

  npy_intp outdim[1] = { 3 };

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "dd", kwlist,
                                  &a, &d))
    goto error;

  sout = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!sout)
    goto error;

  ad_to_v(a, d, PyArray_DATA(sout));

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) sout)));

 error:
  PyArray_XDECREF_ERR((PyArrayObject *) sout);

  return(NULL);
}

static PyObject *lfa_v_to_ad (PyObject *self,
                              PyObject *args,
                              PyObject *kwds) {
  static char *kwlist[] = { "s",
                            "flip",
                            NULL };
  PyObject *sarg;
  unsigned char flip = 0;

  PyObject *sarr = NULL;
  double a, d;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O|b", kwlist,
                                  &sarg, &flip))
    goto error;

  sarr = PyArray_FROM_OTF(sarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!sarr)
    goto error;

  if(PyArray_Size(sarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 's' is too small");
    goto error;
  }

  v_to_ad(PyArray_DATA(sarr), flip, &a, &d);

  Py_DECREF(sarr);

  return(Py_BuildValue("dd", a, d));

 error:
  Py_XDECREF(sarr);

  return(NULL);
}

static PyObject *lfa_v_to_az (PyObject *self,
                              PyObject *args,
                              PyObject *kwds) {
  static char *kwlist[] = { "s",
                            "flip",
                            NULL };
  PyObject *sarg;
  unsigned char flip = 0;

  PyObject *sarr = NULL;
  double a, z;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O|b", kwlist,
                                  &sarg, &flip))
    goto error;

  sarr = PyArray_FROM_OTF(sarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!sarr)
    goto error;

  if(PyArray_Size(sarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 's' is too small");
    goto error;
  }

  v_to_az(PyArray_DATA(sarr), flip, &a, &z);

  Py_DECREF(sarr);

  return(Py_BuildValue("dd", a, z));

 error:
  Py_XDECREF(sarr);

  return(NULL);
}

static PyObject *lfa_v_to_ad_dt (PyObject *self,
                                 PyObject *args,
                                 PyObject *kwds) {
  static char *kwlist[] = { "s",
                            "dsdt",
                            "flip",
                            NULL };
  PyObject *sarg, *dsdtarg;
  unsigned char flip = 0;

  PyObject *sarr = NULL, *dsdtarr = NULL;
  double dadt, dddt;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OO|b", kwlist,
                                  &sarg, &dsdtarg, &flip))
    goto error;

  sarr = PyArray_FROM_OTF(sarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!sarr)
    goto error;

  if(PyArray_Size(sarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 's' is too small");
    goto error;
  }

  dsdtarr = PyArray_FROM_OTF(dsdtarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!dsdtarr)
    goto error;

  if(PyArray_Size(dsdtarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'dsdt' is too small");
    goto error;
  }

  v_to_ad_dt(PyArray_DATA(sarr), PyArray_DATA(dsdtarr), flip,
             &dadt, &dddt);

  Py_DECREF(sarr);
  Py_DECREF(dsdtarr);

  return(Py_BuildValue("dd", dadt, dddt));

 error:
  Py_XDECREF(sarr);
  Py_XDECREF(dsdtarr);

  return(NULL);
}

static PyObject *lfa_v_angle_v (PyObject *self,
                                PyObject *args,
                                PyObject *kwds) {
  static char *kwlist[] = { "u",
                            "v",
                            NULL };
  PyObject *uarg, *varg;

  PyObject *uarr = NULL, *varr = NULL;
  double angle;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OO", kwlist,
                                  &uarg, &varg))
    goto error;

  uarr = PyArray_FROM_OTF(uarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!uarr)
    goto error;

  if(PyArray_Size(uarr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'u' is too small");
    goto error;
  }

  varr = PyArray_FROM_OTF(varg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!varr)
    goto error;

  if(PyArray_Size(varr) < 3) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'v' is too small");
    goto error;
  }

  angle = v_angle_v(PyArray_DATA(uarr), PyArray_DATA(varr));

  Py_DECREF(uarr);
  Py_DECREF(varr);

  return(Py_BuildValue("d", angle));

 error:
  Py_XDECREF(uarr);
  Py_XDECREF(varr);

  return(NULL);
}

static PyObject *lfa_date2mjd (PyObject *self,
                               PyObject *args,
                               PyObject *kwds) {
  static char *kwlist[] = { "yr", "mn", "dy", NULL };
  int yr, mn, dy, mjd;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "iii", kwlist,
                                  &yr, &mn, &dy))
    goto error;

  mjd = date2mjd(yr, mn, dy);

  return(Py_BuildValue("i", mjd));

 error:
  return(NULL);
}

static PyObject *lfa_mjd2date (PyObject *self,
                               PyObject *args,
                               PyObject *kwds) {
  static char *kwlist[] = { "mjd", NULL };
  int mjd, yr, mn, dy;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "i", kwlist,
                                  &mjd))
    goto error;

  mjd2date(mjd, &yr, &mn, &dy);

  return(Py_BuildValue("iii", yr, mn, dy));

 error:
  return(NULL);
}

static PyObject *lfa_mount_ab2rp (PyObject *self,
                                  PyObject *args,
                                  PyObject *kwds) {
  static char *kwlist[] = { "aim",
                            "bore",
                            "snp", "cnp",
                            "flip",
                            "daimdt",
                            NULL };
  PyObject *aimarg, *borearg;
  PyObject *daimdtarg = NULL;
  double snp = 0.0, cnp = 1.0;
  unsigned char flip = 0;

  PyObject *aimarr = NULL;
  PyObject *borearr = NULL;
  PyObject *daimdtarr = NULL;

  double pos[3][3], dposdt[3][3];
  double r, p, drdt, dpdt;

  npy_intp matdim[2] = { 3, 3 };

  PyObject *posout = NULL;
  PyObject *dposdtout = NULL;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OO|ddbO", kwlist,
                                  &aimarg,
                                  &borearg,
                                  &snp, &cnp,
                                  &flip,
                                  &daimdtarg))
    goto error;

  aimarr = PyArray_FROM_OTF(aimarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!aimarr)
    goto error;

  borearr = PyArray_FROM_OTF(borearg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!borearr)
    goto error;

  posout = PyArray_SimpleNew(2, matdim, NPY_DOUBLE);
  if(!posout)
    goto error;

  if(daimdtarg) {
    daimdtarr = PyArray_FROM_OTF(daimdtarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!daimdtarr)
      goto error;

    dposdtout = PyArray_SimpleNew(2, matdim, NPY_DOUBLE);
    if(!dposdtout)
      goto error;
  }

  mount_ab2rp(PyArray_DATA(aimarr),
              daimdtarr ? PyArray_DATA(daimdtarr) : NULL,
              PyArray_DATA(borearr),
              snp, cnp, flip,
              pos, &r, &p,
              dposdt, daimdtarr ? &drdt : NULL, daimdtarr ? &dpdt : NULL);

  memcpy(PyArray_DATA(posout), pos, sizeof(pos));

  Py_DECREF(aimarr);
  Py_DECREF(borearr);

  if(daimdtarr) {
    memcpy(PyArray_DATA(dposdtout), dposdt, sizeof(dposdt));

    Py_DECREF(daimdtarr);

    return(Py_BuildValue("NddNdd",
                         PyArray_Return((PyArrayObject *) posout),
                         r, p,
                         PyArray_Return((PyArrayObject *) dposdtout),
                         drdt, dpdt));
  }
  /* else */

  return(Py_BuildValue("Ndd",
                       PyArray_Return((PyArrayObject *) posout),
                       r, p));

 error:
  Py_XDECREF(aimarr);
  Py_XDECREF(borearr);
  Py_XDECREF(daimdtarr);

  PyArray_XDECREF_ERR((PyArrayObject *) posout);
  PyArray_XDECREF_ERR((PyArrayObject *) dposdtout);

  return(NULL);
}

static PyObject *lfa_mount_pa (PyObject *self,
                               PyObject *args,
                               PyObject *kwds) {
  static char *kwlist[] = { "aimp",
                            "bore",
                            "pos",
                            "daimdt", "dposdt",
                            NULL };
  PyObject *aimparg, *borearg, *posarg;
  PyObject *daimpdtarg = NULL;
  PyObject *dposdtarg = NULL;

  PyObject *aimparr = NULL;
  PyObject *borearr = NULL;
  PyObject *posarr = NULL;
  PyObject *daimpdtarr = NULL;
  PyObject *dposdtarr = NULL;

  double pos[3][3], dposdt[3][3];
  double a, dadt;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OOO|OO", kwlist,
                                  &aimparg,
                                  &borearg,
                                  &posarg,
                                  &daimpdtarg,
                                  &dposdtarg))
    goto error;

  aimparr = PyArray_FROM_OTF(aimparg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!aimparr)
    goto error;

  borearr = PyArray_FROM_OTF(borearg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!borearr)
    goto error;

  posarr = PyArray_FROM_OTF(posarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!posarr)
    goto error;

  memcpy(pos, PyArray_DATA(posarr), sizeof(pos));

  if(daimpdtarg) {
    daimpdtarr = PyArray_FROM_OTF(daimpdtarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!daimpdtarr)
      goto error;

    dposdtarr = PyArray_FROM_OTF(dposdtarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!dposdtarr)
      goto error;

    memcpy(dposdt, PyArray_DATA(dposdtarr), sizeof(dposdt));

    mount_pa(PyArray_DATA(aimparr),
             PyArray_DATA(daimpdtarr),
             PyArray_DATA(borearr),
             pos, dposdt, &a, &dadt);

    Py_DECREF(aimparr);
    Py_DECREF(borearr);
    Py_DECREF(posarr);
    Py_DECREF(daimpdtarr);
    Py_DECREF(dposdtarr);

    return(Py_BuildValue("dd", a, dadt));
  }
  else {
    mount_pa(PyArray_DATA(aimparr),
             NULL,
             PyArray_DATA(borearr),
             pos, dposdt, &a, NULL);

    Py_DECREF(aimparr);
    Py_DECREF(borearr);
    Py_DECREF(posarr);

    return(Py_BuildValue("d", a));
  }

 error:
  Py_XDECREF(aimparr);
  Py_XDECREF(borearr);
  Py_XDECREF(posarr);
  Py_XDECREF(daimpdtarr);
  Py_XDECREF(dposdtarr);

  return(NULL);
}

static PyObject *lfa_mount_rp2pos (PyObject *self,
                                   PyObject *args,
                                   PyObject *kwds) {
  static char *kwlist[] = { "r",
                            "p",
                            "snp",
                            "cnp",
                            NULL };
  double r, p, snp = 0.0, cnp = 1.0;
  double pos[3][3];

  npy_intp matdim[2] = { 3, 3 };
  PyObject *posout = NULL;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "dd|dd", kwlist,
                                  &r, &p, &snp, &cnp))
    goto error;

  posout = PyArray_SimpleNew(2, matdim, NPY_DOUBLE);
  if(!posout)
    goto error;

  mount_rp2pos(r, p, snp, cnp, pos);
  memcpy(PyArray_DATA(posout), pos, sizeof(pos));

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) posout)));

 error:
  PyArray_XDECREF_ERR((PyArrayObject *) posout);

  return(NULL);
}

static PyObject *lfa_dplate (PyObject *self,
                             PyObject *args,
                             PyObject *kwds) {
  static char *kwlist[] = { "comx", "comy",
                            "refx", "refy",
                            "wt",
                            "ncoeff",
                            NULL };
  PyObject *comxarg, *comyarg, *refxarg, *refyarg;
  PyObject *wtarg = NULL;
  int ncoeff = 6;

  PyObject *comxarr = NULL, *comyarr = NULL;
  PyObject *refxarr = NULL, *refyarr = NULL;
  PyObject *wtarr = NULL;
  int npt;

  npy_intp outdim[1] = { 6 };
  PyObject *trout = NULL;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OOOO|Oi", kwlist,
                                  &comxarg, &comyarg,
                                  &refxarg, &refyarg,
                                  &wtarg,
                                  &ncoeff))
    goto error;

  comxarr = PyArray_FROM_OTF(comxarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!comxarr)
    goto error;

  npt = PyArray_Size(comxarr);

  comyarr = PyArray_FROM_OTF(comyarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!comyarr)
    goto error;

  if(PyArray_Size(comyarr) != npt) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'comy' and 'comx' lengths do not match");
    goto error;
  }

  refxarr = PyArray_FROM_OTF(refxarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!refxarr)
    goto error;

  if(PyArray_Size(refxarr) != npt) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'refx' and 'comx' lengths do not match");
    goto error;
  }

  refyarr = PyArray_FROM_OTF(refyarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!refyarr)
    goto error;

  if(PyArray_Size(refyarr) != npt) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'refy' and 'comx' lengths do not match");
    goto error;
  }

  if(wtarg) {
    wtarr = PyArray_FROM_OTF(wtarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!wtarr)
      goto error;
    
    if(PyArray_Size(wtarr) != npt) {
      PyErr_SetString(PyExc_IndexError,
                      "array 'wt' and 'comx' lengths do not match");
      goto error;
    }
  }

  trout = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!trout)
    goto error;

  dplate(PyArray_DATA(comxarr), 0, sizeof(double),
         PyArray_DATA(comyarr), 0, sizeof(double),
         PyArray_DATA(refxarr), 0, sizeof(double),
         PyArray_DATA(refyarr), 0, sizeof(double),
         wtarr ? PyArray_DATA(wtarr) : NULL, 0, sizeof(double),
         npt, ncoeff, PyArray_DATA(trout));

  Py_DECREF(comxarr);
  Py_DECREF(comyarr);
  Py_DECREF(refxarr);
  Py_DECREF(refyarr);
  Py_XDECREF(wtarr);

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) trout)));

 error:
  Py_XDECREF(comxarr);
  Py_XDECREF(comyarr);
  Py_XDECREF(refxarr);
  Py_XDECREF(refyarr);
  Py_XDECREF(wtarr);
  PyArray_XDECREF_ERR((PyArrayObject *) trout);

  return(NULL);
}

static PyObject *lfa_pixovcirc (PyObject *self,
                                PyObject *args,
                                PyObject *kwds) {
  static char *kwlist[] = { "x", "y", "r", NULL };
  PyObject *xarg, *yarg, *rarg;
  PyObject *xarr = NULL, *yarr = NULL, *rarr = NULL;
  double *x, *y, *r;

  PyObject *outarr = NULL;
  double *out;

  int ipt, npt, nr;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OOO", kwlist,
                                  &xarg, &yarg, &rarg))
    goto error;

  xarr = PyArray_FROM_OTF(xarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!xarr)
    goto error;

  npt = PyArray_Size(xarr);

  yarr = PyArray_FROM_OTF(yarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!yarr)
    goto error;

  if(PyArray_Size(yarr) != npt) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'y' and 'x' lengths do not match");
    goto error;
  }

  rarr = PyArray_FROM_OTF(rarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!rarr)
    goto error;

  nr = PyArray_Size(rarr);

  if(nr != 1 && nr != npt) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'r' and 'x' lengths do not match");
    goto error;
  }

  outarr = PyArray_SimpleNew(PyArray_NDIM(xarr),
                             PyArray_DIMS(xarr),
                             NPY_DOUBLE);
  if(!outarr)
    goto error;

  x = PyArray_DATA(xarr);
  y = PyArray_DATA(yarr);
  r = PyArray_DATA(rarr);
  out = PyArray_DATA(outarr);

  for(ipt = 0; ipt < npt; ipt++)
    out[ipt] = pixovcirc(x[ipt], y[ipt], nr == 1 ? r[0] : r[ipt]);

  Py_DECREF(xarr);
  Py_DECREF(yarr);
  Py_DECREF(rarr);

  return(PyArray_Return((PyArrayObject *) outarr));

 error:
  Py_XDECREF(xarr);
  Py_XDECREF(yarr);
  Py_XDECREF(rarr);
  PyArray_XDECREF_ERR((PyArrayObject *) outarr);

  return(NULL);
}

static PyObject *lfa_skylevel_image (struct lfa_ap_object *self,
                                     PyObject *args,
                                     PyObject *kwds) {
  static char *kwlist[] = { "map",
                            "mask",
                            "clip_low",
                            "clip_high",
                            NULL };
  PyObject *maparg;
  PyObject *maskarg = NULL;
  float clip_low = -FLT_MAX;
  float clip_high = 3.0;
  float skylev, skynoise;

  PyObject *maparr = NULL;
  PyObject *maskarr = NULL;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O|Off", kwlist,
                                  &maparg,
                                  &maskarg,
                                  &clip_low,
                                  &clip_high))
    return(NULL);

  maparr = PyArray_FROM_OTF(maparg, NPY_FLOAT, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!maparr)
    goto error;

  if(maskarg) {
    maskarr = PyArray_FROM_OTF(maskarg, NPY_UINT8, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!maskarr)
      goto error;

    if(PyArray_Size(maskarr) < PyArray_Size(maparr)) {
      PyErr_SetString(PyExc_IndexError, "mask too small");
      goto error;
    }
  }

  skylevel_image(PyArray_DATA(maparr),
                 maskarr ? PyArray_DATA(maskarr) : NULL,
                 PyArray_Size(maparr),
                 clip_low, clip_high,
                 &skylev, &skynoise);

  Py_DECREF(maparr);
  Py_XDECREF(maskarr);

  return(Py_BuildValue("ff", skylev, skynoise));

 error:
  Py_XDECREF(maparr);
  Py_XDECREF(maskarr);

  return(NULL);
}

static PyMethodDef lfa_methods[] = {
  { "source_star_vec", (PyCFunction) lfa_source_star_vec,
    METH_VARARGS | METH_KEYWORDS,
    "src = source_star_vec(n, dndt=None, pr=0, epoch=2000.0)" },
  { "source_elem", (PyCFunction) lfa_source_elem,
    METH_VARARGS | METH_KEYWORDS,
    "src = source_elem(eltype, epoch, incl, anode, longperi, aq, ecc, lm, nn)" },
  { "source_mpc", (PyCFunction) lfa_source_mpc,
    METH_VARARGS | METH_KEYWORDS,
    "src = source_mpc(line)" },
  { "source_planet", (PyCFunction) lfa_source_planet,
    METH_VARARGS | METH_KEYWORDS,
    "src = source_planet(n)" },
  { "base60_to_10", (PyCFunction) lfa_base60_to_10,
    METH_VARARGS | METH_KEYWORDS,
    "outangle = base60_to_10(inangle, sep, inunit, outunit)" },
  { "base10_to_60", (PyCFunction) lfa_base10_to_60,
    METH_VARARGS | METH_KEYWORDS,
    "outangle = base10_to_60(inangle, inunit, sep, sign, dp, outunit)" },
  { "vec2tp", (PyCFunction) lfa_vec2tp,
    METH_VARARGS | METH_KEYWORDS,
    "x, y = vec2tp(s, tp)" },
  { "tp2vec", (PyCFunction) lfa_tp2vec,
    METH_VARARGS | METH_KEYWORDS,
    "s = tp2vec(x, y, tp)" },
  { "v_airmass", (PyCFunction) lfa_v_airmass,
    METH_VARARGS | METH_KEYWORDS,
    "airmass = v_airmass(s)" },
  { "v_parallactic", (PyCFunction) lfa_v_parallactic,
    METH_VARARGS | METH_KEYWORDS,
    "pa = v_parallactic(sinphi, cosphi, s)" },
  { "ad_to_v", (PyCFunction) lfa_ad_to_v,
    METH_VARARGS | METH_KEYWORDS,
    "s = ad_to_v(a, d)" },
  { "v_to_ad", (PyCFunction) lfa_v_to_ad,
    METH_VARARGS | METH_KEYWORDS,
    "a, d = v_to_ad(s, flip=0)" },
  { "v_to_az", (PyCFunction) lfa_v_to_az,
    METH_VARARGS | METH_KEYWORDS,
    "a, z = v_to_az(s, flip=0)" },
  { "v_to_ad_dt", (PyCFunction) lfa_v_to_ad_dt,
    METH_VARARGS | METH_KEYWORDS,
    "dadt, dddt = v_to_ad_dt(s, dsdt, flip=0)" },
  { "v_angle_v", (PyCFunction) lfa_v_angle_v,
    METH_VARARGS | METH_KEYWORDS,
    "angle = v_angle_v(u, v)" },
  { "date2mjd", (PyCFunction) lfa_date2mjd,
    METH_VARARGS | METH_KEYWORDS,
    "mjd = date2mjd(yr, mn, dy)" },
  { "mjd2date", (PyCFunction) lfa_mjd2date,
    METH_VARARGS | METH_KEYWORDS,
    "yr, mn, dy = mjd2date(mjd)" },
  { "mount_ab2rp", (PyCFunction) lfa_mount_ab2rp,
    METH_VARARGS | METH_KEYWORDS,
    "pos, r, p = mount_rp2pos(aim, bore, snp=0, cnp=1, flip=0, daimdt=None)" },
  { "mount_pa", (PyCFunction) lfa_mount_pa,
    METH_VARARGS | METH_KEYWORDS,
    "a = mount_pa(aimp, bore, pos)" },
  { "mount_rp2pos", (PyCFunction) lfa_mount_rp2pos,
    METH_VARARGS | METH_KEYWORDS,
    "pos = mount_rp2pos(r, p, snp=0, cnp=1)" },
  { "dplate", (PyCFunction) lfa_dplate,
    METH_VARARGS | METH_KEYWORDS,
    "tr = dplate(comx, comy, refx, refy, wt=None, ncoeff=6)" },
  { "pixovcirc", (PyCFunction) lfa_pixovcirc,
    METH_VARARGS | METH_KEYWORDS,
    "fract = pixovcirc(x, y, r)" },
  { "skylevel_image", (PyCFunction) lfa_skylevel_image,
    METH_VARARGS | METH_KEYWORDS,
    "skylev, skynoise = skylevel_image(map, mask=None, clip_low=-FLT_MAX, clip_high=3)" },
  { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC initlfa (void) {
  PyObject *m, *o;

#define MAKECONST(m) { #m, m }

  struct {
    char *name;
    int value;
  } iconst[] = {
    MAKECONST(TR_MOTION),
    MAKECONST(TR_PLX),
    MAKECONST(TR_DEFL),
    MAKECONST(TR_ANNAB),
    MAKECONST(TR_TOPO),
    MAKECONST(TR_LAT),
    MAKECONST(TR_REFRO),
    MAKECONST(TR_TO_AST),
    MAKECONST(TR_TO_GCRS),
    MAKECONST(TR_TO_TOPO_HD),
    MAKECONST(TR_TO_OBS_HD),
    MAKECONST(TR_TO_TOPO_AZ),
    MAKECONST(TR_TO_OBS_AZ),

    MAKECONST(JPLEPH_MERCURY),
    MAKECONST(JPLEPH_VENUS),
    MAKECONST(JPLEPH_EMB),
    MAKECONST(JPLEPH_MARS),
    MAKECONST(JPLEPH_JUPITER),
    MAKECONST(JPLEPH_SATURN),
    MAKECONST(JPLEPH_URANUS),
    MAKECONST(JPLEPH_NEPTUNE),
    MAKECONST(JPLEPH_PLUTO),
    MAKECONST(JPLEPH_MOON),
    MAKECONST(JPLEPH_SUN),
    MAKECONST(JPLEPH_NUTATION),
    MAKECONST(JPLEPH_LIBRATION),
    MAKECONST(JPLEPH_EULER),
    MAKECONST(JPLEPH_TEI),
    MAKECONST(TIMEEPH_TEI),
    MAKECONST(TIMEEPH_TEV),
    MAKECONST(SOURCE_STAR),
    MAKECONST(SOURCE_ELEM),

    MAKECONST(OBSERVER_UPDATE_IERS),
    MAKECONST(OBSERVER_UPDATE_NUT),
    MAKECONST(OBSERVER_UPDATE_PFB),
    MAKECONST(OBSERVER_UPDATE_ERA),
    MAKECONST(OBSERVER_UPDATE_TDB),
    MAKECONST(OBSERVER_UPDATE_SOLSYS),
    MAKECONST(OBSERVER_UPDATE_ALL),

    MAKECONST(UNIT_RAD),
    MAKECONST(UNIT_HR),
    MAKECONST(UNIT_M),
    MAKECONST(UNIT_S),
    MAKECONST(UNIT_MS),
    MAKECONST(UNIT_DEG),
    MAKECONST(UNIT_AM),
    MAKECONST(UNIT_AS),
    MAKECONST(UNIT_MAS)
  };

  struct {
    char *name;
    double value;
  } dconst[] = {
    MAKECONST(GMSUN),
    MAKECONST(AU),
    MAKECONST(LIGHT),
    MAKECONST(GMEARTH),
    MAKECONST(AEARTH),
    MAKECONST(FEARTH),
    MAKECONST(DAY),
    MAKECONST(DTT),
    MAKECONST(ZMJD),
    MAKECONST(JYR),
    MAKECONST(J2K),
    MAKECONST(JUNIX),
    MAKECONST(JTCB),
    MAKECONST(LB),
    MAKECONST(LG),
    MAKECONST(ZTDB),
    MAKECONST(TWOPI),
    MAKECONST(AS_PER_REV),
    MAKECONST(ERA2K),
    MAKECONST(ERADAY),
    MAKECONST(EOMEGA),
    MAKECONST(RSUN),
    MAKECONST(RJUP),
    MAKECONST(REARTH),

    MAKECONST(MAS_TO_RAD),
    MAKECONST(DEG_TO_RAD),
    MAKECONST(AM_TO_RAD),
    MAKECONST(AS_TO_RAD),
    MAKECONST(HR_TO_RAD),
    MAKECONST(SEC_TO_RAD),
    MAKECONST(RAD_TO_DEG),
    MAKECONST(RAD_TO_AM),
    MAKECONST(RAD_TO_AS),
    MAKECONST(RAD_TO_HR),
    MAKECONST(DEG_TO_AS)
  };

  int c, nc;

  npy_intp matdim[2] = { 3, 3 };

  /* Init module */
  m = Py_InitModule("lfa", lfa_methods);
  if(!m)
    return;

  /* Import numpy */
  import_array();

  /* Create constants */
  nc = sizeof(iconst) / sizeof(iconst[0]);

  for(c = 0; c < nc; c++) {
    o = PyInt_FromLong(iconst[c].value);
    if(o)
      PyModule_AddObject(m, iconst[c].name, o);
  }

  nc = sizeof(dconst) / sizeof(dconst[0]);

  for(c = 0; c < nc; c++) {
    o = PyFloat_FromDouble(dconst[c].value);
    if(o)
      PyModule_AddObject(m, dconst[c].name, o);
  }

  /* Constant matrices */
#define MAKEMATRIX(mat)                                         \
  o = PyArray_SimpleNewFromData(2, matdim, NPY_DOUBLE, mat);    \
  if(o)                                                         \
    PyModule_AddObject(m, #mat, o)

  MAKEMATRIX(gcrs2ecl);
  MAKEMATRIX(ecl2gcrs);
  MAKEMATRIX(fk52ecl);
  MAKEMATRIX(ecl2fk5);
  MAKEMATRIX(eq2gal);
  MAKEMATRIX(gal2eq);

  /* Create types */
  lfa_source_type.tp_new = PyType_GenericNew;
  if(PyType_Ready(&lfa_source_type) < 0)
    return;

  Py_INCREF(&lfa_source_type);
  PyModule_AddObject(m, "source", (PyObject *) &lfa_source_type);

  if(PyType_Ready(&lfa_observer_type) < 0)
    return;

  Py_INCREF(&lfa_observer_type);
  PyModule_AddObject(m, "observer", (PyObject *) &lfa_observer_type);

  if(PyType_Ready(&lfa_ap_type) < 0)
    return;

  Py_INCREF(&lfa_ap_type);
  PyModule_AddObject(m, "ap", (PyObject *) &lfa_ap_type);
}
