#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>

#include <stdlib.h>
#include <inttypes.h>
#include <float.h>

#include "ap.h"
#include "cvtunit.h"
#include "specfind.h"
#include "lfa.h"

/* This macro, which was recommended by all of the documentation back
 * when I wrote this module, was deprecated in Numpy API 1.14 and
 * removed in 1.24.  There seems to be no good reason for removing it
 * as far as I can tell, especially when the recommended replacement
 * is two lines rather than one and uses a function that was brand new
 * in 1.14 so isn't backwards compatible.  We therefore just restore
 * the macro using the currently recommended method so we can retain
 * compatibility with pre-1.14 numpy API. */
#ifndef PyArray_XDECREF_ERR
#define PyArray_XDECREF_ERR(o)                                  \
  PyArray_DiscardWritebackIfCopy((PyArrayObject *) (o));        \
  Py_XDECREF(o)
#endif

/* If we're compiling without the pre-1.7 API, restore these macros */
#ifndef NPY_FORCECAST
#define NPY_FORCECAST NPY_ARRAY_FORCECAST
#endif

#ifndef NPY_IN_ARRAY
#define NPY_IN_ARRAY NPY_ARRAY_IN_ARRAY
#endif

/* Macros to help me not get confused about y,x in numpy 2D arrays */
#define DIM_Y  0
#define DIM_X  1

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
  PyArray_Descr *dtype;
};

struct lfa_fsgp_object {
  PyObject_HEAD
  struct fsgp_fac fac;
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
  PyVarObject_HEAD_INIT(NULL, 0)
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

  Py_TYPE(self)->tp_free((PyObject *) self);
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

  Py_INCREF(Py_None);
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
  PyVarObject_HEAD_INIT(NULL, 0)
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
    Py_DECREF(self->dtype);
    self->allocated = 0;
  }

  Py_TYPE(self)->tp_free((PyObject *) self);
}

/* C preprocessor stringification of a macro requires an additional
   layer of expansion, so we need to define two macros to achieve it. */
#undef mstrify
#undef strify

#define mstrify(x) strify(x)
#define strify(x) #x

static PyArray_Descr *lfa_ap_build_descr (void) {
  PyObject *names = (PyObject *) NULL;
  PyObject *formats = (PyObject *) NULL;
  PyObject *offsets = (PyObject *) NULL;
  PyObject *d = (PyObject *) NULL;
  PyArray_Descr *dtype = (PyArray_Descr *) NULL;
  PyObject *o;

  struct {
    char *name;
    char *format;
    size_t offset;
  } members[] = {
    { "x", "f8", offsetof(struct ap_source, x) },
    { "y", "f8", offsetof(struct ap_source, y) },
    { "sxx", "f8", offsetof(struct ap_source, sxx) },
    { "syy", "f8", offsetof(struct ap_source, syy) },
    { "sxy", "f8", offsetof(struct ap_source, sxy) },
    { "ziso", "f8", offsetof(struct ap_source, ziso) },
    { "zpk", "f8", offsetof(struct ap_source, zpk) },
    { "hfd", "f8", offsetof(struct ap_source, hfd) },
    { "areal", mstrify(AP_NAREAL) "i4", offsetof(struct ap_source, areal) },
    { "bb", "4i4", offsetof(struct ap_source, areal) },
    { "flags", "B", offsetof(struct ap_source, flags) },
  };

  int imemb, nmemb;

  /* Create lists */
  nmemb = sizeof(members) / sizeof(members[0]);

  names = PyList_New(nmemb);
  formats = PyList_New(nmemb);
  offsets = PyList_New(nmemb);
  if(!names || !formats || !offsets)
    goto error;

  /* Populate */
  for(imemb = 0; imemb < nmemb; imemb++) {
#if PY_MAJOR_VERSION >= 3
    o = PyUnicode_FromString(members[imemb].name);
#else
    o = PyString_FromString(members[imemb].name);
#endif
    if(!o)
      goto error;

    PyList_SetItem(names, imemb, o);  /* list takes over our ref to o */

#if PY_MAJOR_VERSION >= 3
    o = PyUnicode_FromString(members[imemb].format);
#else
    o = PyString_FromString(members[imemb].format);
#endif
    if(!o)
      goto error;

    PyList_SetItem(formats, imemb, o);

#if PY_MAJOR_VERSION >= 3
    o = PyLong_FromSize_t(members[imemb].offset);
#else
    o = PyInt_FromSize_t(members[imemb].offset);
#endif
    if(!o)
      goto error;

    PyList_SetItem(offsets, imemb, o);
  }

  /* Create dictionary */
  d = PyDict_New();
  if(!d)
    goto error;

  PyDict_SetItemString(d, "names", names);  /* does not take over our ref */
  Py_DECREF(names);                         /* so we have to do this */

  PyDict_SetItemString(d, "formats", formats);
  Py_DECREF(formats);

  PyDict_SetItemString(d, "offsets", offsets);
  Py_DECREF(offsets);

#if PY_MAJOR_VERSION >= 3
  o = PyLong_FromSize_t(sizeof(struct ap_source));
#else
  o = PyInt_FromSize_t(sizeof(struct ap_source));
#endif
  if(!o)
    goto error;

  PyDict_SetItemString(d, "itemsize", o);
  Py_DECREF(o);

  /* Convert to descriptor */
  if(PyArray_DescrConverter(d, &dtype) != NPY_SUCCEED)
    goto error;

  Py_DECREF(d);

  return(dtype);

 error:
  Py_XDECREF(names);
  Py_XDECREF(formats);
  Py_XDECREF(offsets);
  Py_XDECREF(d);

  return(NULL);
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

  self->dtype = lfa_ap_build_descr();
  if(!self->dtype) {
    ap_free(&(self->ap));
    return(-1);
  }

  self->allocated = 1;

  return(0);
}

/* In order to avoid leaking source lists when we wrap them in numpy
   arrays, we need to create an object to "own" the source list. */

struct lfa_ap_source_list_object {
  PyObject_HEAD
  struct ap_source *list;
  int nlist;
};

static PyObject *lfa_ap_source_list_new (PyTypeObject *type,
                                         PyObject *args,
                                         PyObject *kwds) {
  struct lfa_ap_source_list_object *self = NULL;

  /* Allocate */
  self = (struct lfa_ap_source_list_object *) type->tp_alloc(type, 0);
  if(self) {
    self->list = NULL;
    self->nlist = 0;
  }

  return((PyObject *) self);
}

static void lfa_ap_source_list_dealloc (struct lfa_ap_source_list_object *self) {
  if(self->list) {
    free(self->list);
    self->list = NULL;
    self->nlist = 0;
  }

  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyTypeObject lfa_ap_source_list_type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "lfa.ap.source_list",                      /* tp_name */
  sizeof(struct lfa_ap_source_list_object),  /* tp_basicsize */
  0,                                         /* tp_itemsize */
  (destructor) lfa_ap_source_list_dealloc,   /* tp_dealloc */
  0,                                         /* tp_print */
  0,                                         /* tp_getattr */
  0,                                         /* tp_setattr */
  0,                                         /* tp_compare */
  0,                                         /* tp_repr */
  0,                                         /* tp_as_number */
  0,                                         /* tp_as_sequence */
  0,                                         /* tp_as_mapping */
  0,                                         /* tp_hash */
  0,                                         /* tp_call */
  0,                                         /* tp_str */
  0,                                         /* tp_getattro */
  0,                                         /* tp_setattro */
  0,                                         /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,  /* tp_flags */
  "ap source list",                          /* tp_doc */
  0,                                         /* tp_traverse */
  0,                                         /* tp_clear */
  0,                                         /* tp_richcompare */
  0,                                         /* tp_weaklistoffset */
  0,                                         /* tp_iter */
  0,                                         /* tp_iternext */
  0,                                         /* tp_methods */
  0,                                         /* tp_members */
  0,                                         /* tp_getset */
  0,                                         /* tp_base */
  0,                                         /* tp_dict */
  0,                                         /* tp_descr_get */
  0,                                         /* tp_descr_set */
  0,                                         /* tp_dictoffset */
  0,                                         /* tp_init */
  0,                                         /* tp_alloc */
  lfa_ap_source_list_new                     /* tp_new */
};

static int lfa_ap_image_output (struct ap *ap, struct ap_source *obj,
                                void *data) {
  struct lfa_ap_source_list_object *h
    = (struct lfa_ap_source_list_object *) data;

  h->list = realloc(h->list, (h->nlist+1) * sizeof(struct ap_source));
  if(!h->list)
    return(-1);

  memcpy(&(h->list[h->nlist]), obj, sizeof(struct ap_source));
  h->nlist++;

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

  struct lfa_ap_source_list_object *slist = NULL;
  npy_intp outdim[1] = { 0 };
  PyObject *result = NULL;

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

  if(filtmaparg && filtmaparg != Py_None) {
    filtmaparr = PyArray_FROM_OTF(filtmaparg, NPY_FLOAT, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!filtmaparr)
      goto error;

    if(PyArray_Size(filtmaparr) < self->ap.nx * self->ap.ny) {
      PyErr_SetString(PyExc_IndexError, "filtmap too small");
      goto error;
    }
  }

  if(maskarg && maskarg != Py_None) {
    maskarr = PyArray_FROM_OTF(maskarg, NPY_UINT8, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!maskarr)
      goto error;

    if(PyArray_Size(maskarr) < self->ap.nx * self->ap.ny) {
      PyErr_SetString(PyExc_IndexError, "mask too small");
      goto error;
    }
  }

  /* Create source list object */
  slist = (struct lfa_ap_source_list_object *)
    lfa_ap_source_list_type.tp_alloc(&lfa_ap_source_list_type, 0);
  if(!slist)
    goto error;

  slist->list = NULL;
  slist->nlist = 0;

  /* Hook it into output */
  self->ap.output = lfa_ap_image_output;
  self->ap.output_user_data = (void *) slist;

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

  /* Create result array.  This steals a reference to "dtype",
     so we need to get one to hand over.  The memory buffer for
     the source list is owned by "slist" */
  outdim[0] = slist->nlist;

  Py_INCREF(self->dtype);
  result = PyArray_NewFromDescr(&PyArray_Type,
                                self->dtype,
                                1, outdim, NULL,
                                slist->list, 0,
                                NULL);
  if(!result) {
    Py_DECREF(self->dtype);
    goto error;
  }

  /* Tell numpy who owns the memory.  This steals our reference to the
     source list, so we don't need to decref it. */
#if defined(NPY_API_VERSION) && defined(NPY_1_7_API_VERSION) && NPY_API_VERSION >= NPY_1_7_API_VERSION
  PyArray_SetBaseObject((PyArrayObject *) result, (PyObject *) slist);
#else
  PyArray_BASE((PyArrayObject *) result) = (PyObject *) slist;
#endif

  /* Finally we can return the result */
  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) result)));

 error:
  Py_XDECREF(maparr);
  Py_XDECREF(filtmaparr);
  Py_XDECREF(maskarr);
  Py_XDECREF((PyObject *) slist);
  PyArray_XDECREF_ERR((PyArrayObject *) result);

  return(NULL);
}

static PyObject *lfa_ap_ellipse (struct lfa_ap_object *self,
                                 PyObject *args,
                                 PyObject *kwds) {
  static char *kwlist[] = { "source",
                            NULL };

  PyObject *srcarg;
  PyObject *srcarr = NULL;
  struct ap_source *srclist;

  PyObject *gfwhmout = NULL;
  PyObject *ellout = NULL;
  PyObject *spaout = NULL;
  PyObject *cpaout = NULL;
  double *gfwhm, *ell, *spa, *cpa;

  int i, n;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O", kwlist,
                                  &srcarg))
    return(NULL);

  /* Extract array.  This steals a reference to "dtype",
     so we need to get one to hand over. */
  Py_INCREF(self->dtype);
  srcarr = PyArray_FromAny(srcarg, self->dtype, 0, 0,
                           NPY_IN_ARRAY | NPY_FORCECAST, NULL);
  if(!srcarr) {
    Py_DECREF(self->dtype);
    goto error;
  }

  /* Create outputs */
  gfwhmout = PyArray_SimpleNew(PyArray_NDIM(srcarr),
                               PyArray_DIMS(srcarr),
                               NPY_DOUBLE);
  if(!gfwhmout)
    goto error;

  ellout = PyArray_SimpleNew(PyArray_NDIM(srcarr),
                             PyArray_DIMS(srcarr),
                             NPY_DOUBLE);
  if(!ellout)
    goto error;

  spaout = PyArray_SimpleNew(PyArray_NDIM(srcarr),
                               PyArray_DIMS(srcarr),
                               NPY_DOUBLE);
  if(!spaout)
    goto error;

  cpaout = PyArray_SimpleNew(PyArray_NDIM(srcarr),
                               PyArray_DIMS(srcarr),
                               NPY_DOUBLE);
  if(!cpaout)
    goto error;

  n = PyArray_SIZE(srcarr);
  srclist = PyArray_DATA(srcarr);

  gfwhm = PyArray_DATA(gfwhmout);
  ell = PyArray_DATA(ellout);
  spa = PyArray_DATA(spaout);
  cpa = PyArray_DATA(cpaout);

  for(i = 0; i < n; i++)
    ap_ellipse(&(srclist[i]),
               &(gfwhm[i]), &(ell[i]), &(spa[i]), &(cpa[i]));

  Py_DECREF(srcarr);

  return(Py_BuildValue("NNNN",
                       PyArray_Return((PyArrayObject *) gfwhmout),
                       PyArray_Return((PyArrayObject *) ellout),
                       PyArray_Return((PyArrayObject *) spaout),
                       PyArray_Return((PyArrayObject *) cpaout)));

 error:
  Py_XDECREF(srcarr);
  PyArray_XDECREF_ERR((PyArrayObject *) gfwhmout);
  PyArray_XDECREF_ERR((PyArrayObject *) ellout);
  PyArray_XDECREF_ERR((PyArrayObject *) spaout);
  PyArray_XDECREF_ERR((PyArrayObject *) cpaout);

  return(NULL);
}

static PyMemberDef lfa_ap_members[] = {
  { NULL, 0, 0, 0, NULL }
};

static PyMethodDef lfa_ap_methods[] = {
  { "image", (PyCFunction) lfa_ap_image,
    METH_VARARGS | METH_KEYWORDS,
    "sources = image(map, minpix, sky, thresh)" },
  { "ellipse", (PyCFunction) lfa_ap_ellipse,
    METH_VARARGS | METH_KEYWORDS,
    "gfwhm, ell, spa, cpa = ellipse(source)" },
  { NULL, NULL, 0, NULL }
};

static PyTypeObject lfa_ap_type = {
  PyVarObject_HEAD_INIT(NULL, 0)
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

static PyObject *lfa_ranorm (PyObject *self,
                             PyObject *args,
                             PyObject *kwds) {
  static char *kwlist[] = { "a", NULL };
  PyObject *aarg;
  PyObject *aarr = NULL;
  double *a;

  PyObject *outarr = NULL;
  double *out;

  int ipt, npt;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O", kwlist,
                                  &aarg))
    goto error;

  aarr = PyArray_FROM_OTF(aarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!aarr)
    goto error;

  npt = PyArray_Size(aarr);

  outarr = PyArray_SimpleNew(PyArray_NDIM(aarr),
                             PyArray_DIMS(aarr),
                             NPY_DOUBLE);
  if(!outarr)
    goto error;

  a = PyArray_DATA(aarr);
  out = PyArray_DATA(outarr);

  for(ipt = 0; ipt < npt; ipt++)
    out[ipt] = ranorm(a[ipt]);

  Py_DECREF(aarr);

  return(PyArray_Return((PyArrayObject *) outarr));

 error:
  Py_XDECREF(aarr);
  PyArray_XDECREF_ERR((PyArrayObject *) outarr);

  return(NULL);
}

static PyObject *lfa_range (PyObject *self,
                            PyObject *args,
                            PyObject *kwds) {
  static char *kwlist[] = { "a", NULL };
  PyObject *aarg;
  PyObject *aarr = NULL;
  double *a;

  PyObject *outarr = NULL;
  double *out;

  int ipt, npt;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O", kwlist,
                                  &aarg))
    goto error;

  aarr = PyArray_FROM_OTF(aarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!aarr)
    goto error;

  npt = PyArray_Size(aarr);

  outarr = PyArray_SimpleNew(PyArray_NDIM(aarr),
                             PyArray_DIMS(aarr),
                             NPY_DOUBLE);
  if(!outarr)
    goto error;

  a = PyArray_DATA(aarr);
  out = PyArray_DATA(outarr);

  for(ipt = 0; ipt < npt; ipt++)
    out[ipt] = range(a[ipt]);

  Py_DECREF(aarr);

  return(PyArray_Return((PyArrayObject *) outarr));

 error:
  Py_XDECREF(aarr);
  PyArray_XDECREF_ERR((PyArrayObject *) outarr);

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

static PyObject *lfa_geoc (PyObject *self,
                           PyObject *args,
                           PyObject *kwds) {
  static char *kwlist[] = { "sinphi",
                            "cosphi",
                            "height",
                            NULL };
  double sinphi, cosphi, height;
  double u, z;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "ddd", kwlist,
                                  &sinphi,
                                  &cosphi,
                                  &height))
    goto error;

  geoc(sinphi, cosphi, height, &u, &z);

  return(Py_BuildValue("dd", u, z));

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

  if(daimdtarg && daimdtarg != Py_None) {
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

  if(daimpdtarg && daimpdtarg != Py_None) {
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

static PyObject *lfa_pfb_matrix (PyObject *self,
                                 PyObject *args,
                                 PyObject *kwds) {
  static char *kwlist[] = { "tt",
                            NULL };
  double tt, jctk;
  double ang[NPNANG];
  double mat[3][3];

  npy_intp matdim[2] = { 3, 3 };
  PyObject *pfbout = NULL;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "d", kwlist,
                                  &tt))
    goto error;

  /* Create output array */
  pfbout = PyArray_SimpleNew(2, matdim, NPY_DOUBLE);
  if(!pfbout)
    goto error;

  /* Time argument: TT Julian centuries since 2000.0 */
  jctk = (tt-J2K) / (100*JYR);

  /* Get Fukushima-Williams angles */
  pfb06ang(jctk, ang);

  /* Form matrix */
  makepnm(ang, mat);

  memcpy(PyArray_DATA(pfbout), mat, sizeof(mat));

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) pfbout)));

 error:
  PyArray_XDECREF_ERR((PyArrayObject *) pfbout);

  return(NULL);
}

static PyObject *lfa_refract_const (PyObject *self,
                                    PyObject *args,
                                    PyObject *kwds) {
  static char *kwlist[] = { "temperat",
                            "humidity",
                            "pressure",
                            "wavelength",
                            "height",
                            NULL };
  double temperat, humidity, pressure, wavelength, height;

  PyObject *out = NULL;

  npy_intp outdim[1] = { NREFCO };

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "ddddd", kwlist,
                                  &temperat,
                                  &humidity,
                                  &pressure,
                                  &wavelength,
                                  &height))
    return(NULL);

  out = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!out)
    goto error;

  refract_const(temperat, humidity, pressure, wavelength,
                height, PyArray_DATA(out));

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) out)));

 error:
  PyArray_XDECREF_ERR((PyArrayObject *) out);

  return(NULL);
}

static PyObject *lfa_refract_corr (PyObject *self,
                                   PyObject *args,
                                   PyObject *kwds) {
  static char *kwlist[] = { "refco",
                            "tanz",
                            NULL };
  PyObject *rarg;

  PyObject *rarr = NULL;

  double tanz, refr, deriv;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "Od", kwlist,
                                  &rarg, &tanz))
    goto error;

  rarr = PyArray_FROM_OTF(rarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!rarr)
    goto error;

  if(PyArray_Size(rarr) < NREFCO) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'refco' is too small");
    goto error;
  }

  refract_corr(PyArray_DATA(rarr), tanz, &refr, &deriv);

  Py_DECREF(rarr);

  return(Py_BuildValue("dd", refr, deriv));

 error:
  Py_XDECREF(rarr);

  return(NULL);
}

static PyObject *lfa_refract_vec (PyObject *self,
                                  PyObject *args,
                                  PyObject *kwds) {
  static char *kwlist[] = { "refco",
                            "vi",
                            "unref",
                            "dvidt",
                            NULL };
  PyObject *rarg, *varg, *dvarg = NULL;
  int unref = 0;

  PyObject *rarr = NULL, *varr = NULL, *dvarr = NULL;

  PyObject *out = NULL, *doutdt = NULL;

  npy_intp outdim[1] = { 3 };

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OO|iO", kwlist,
                                  &rarg, &varg, &unref, &dvarg))
    goto error;

  rarr = PyArray_FROM_OTF(rarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!rarr)
    goto error;

  if(PyArray_Size(rarr) < NREFCO) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'refco' is too small");
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

  if(dvarg && dvarg != Py_None) {
    dvarr = PyArray_FROM_OTF(dvarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!dvarr)
      goto error;

    if(PyArray_Size(dvarr) < 3) {
      PyErr_SetString(PyExc_IndexError,
                      "array 'dvdt' is too small");
      goto error;
    }
  }

  out = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
  if(!out)
    goto error;

  if(dvarr) {
    doutdt = PyArray_SimpleNew(1, outdim, NPY_DOUBLE);
    if(!doutdt)
      goto error;
  }

  refract_vec(PyArray_DATA(rarr), unref,
              PyArray_DATA(varr), PyArray_DATA(out),
              dvarr ? PyArray_DATA(dvarr) : NULL,
              dvarr ? PyArray_DATA(doutdt) : NULL);

  Py_DECREF(rarr);
  Py_DECREF(varr);
  if(dvarr)
    Py_DECREF(dvarr);

  if(dvarr)
    return(Py_BuildValue("NN",
                         PyArray_Return((PyArrayObject *) out),
                         PyArray_Return((PyArrayObject *) doutdt)));
  else
    return(Py_BuildValue("N",
                         PyArray_Return((PyArrayObject *) out)));

 error:
  Py_XDECREF(rarr);
  Py_XDECREF(varr);
  Py_XDECREF(dvarr);
  PyArray_XDECREF_ERR((PyArrayObject *) out);
  PyArray_XDECREF_ERR((PyArrayObject *) doutdt);

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

  if(wtarg && wtarg != Py_None) {
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

static PyObject *lfa_kepler (PyObject *self,
                             PyObject *args,
                             PyObject *kwds) {
  static char *kwlist[] = { "ma", "ecc", NULL };
  PyObject *maarg;
  PyObject *maarr = NULL;
  double *ma;
  double ecc;

  PyObject *eaarr = NULL;
  double *ea;

  int ipt, npt;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "Od", kwlist,
                                  &maarg, &ecc))
    goto error;

  maarr = PyArray_FROM_OTF(maarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!maarr)
    goto error;

  npt = PyArray_Size(maarr);

  eaarr = PyArray_SimpleNew(PyArray_NDIM(maarr),
                            PyArray_DIMS(maarr),
                            NPY_DOUBLE);
  if(!eaarr)
    goto error;

  ma = PyArray_DATA(maarr);
  ea = PyArray_DATA(eaarr);

  for(ipt = 0; ipt < npt; ipt++)
    ea[ipt] = kepler(ma[ipt], ecc);

  Py_DECREF(maarr);

  return(PyArray_Return((PyArrayObject *) eaarr));

 error:
  Py_XDECREF(maarr);
  PyArray_XDECREF_ERR((PyArrayObject *) eaarr);

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

static PyObject *lfa_backremove (PyObject *self,
                                 PyObject *args,
                                 PyObject *kwds) {
  static char *kwlist[] = { "map",
                            "mask",
                            "nbsize",
                            NULL };
  PyObject *maparg;
  PyObject *maskarg = NULL;
  int nbsize = 64;

  PyObject *maparr = NULL;
  PyObject *maskarr = NULL;
  PyObject *outarr = NULL;

  int ndim;
  npy_intp *dims;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O|Oi", kwlist,
                                  &maparg,
                                  &maskarg,
                                  &nbsize))
    return(NULL);

  maparr = PyArray_FROM_OTF(maparg, NPY_FLOAT, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!maparr)
    goto error;

  ndim = PyArray_NDIM(maparr);
  dims = PyArray_DIMS(maparr);

  if(ndim != 2) {
    PyErr_SetString(PyExc_IndexError, "expected 2-D map");
    goto error;
  }

  if(maskarg && maskarg != Py_None) {
    maskarr = PyArray_FROM_OTF(maskarg, NPY_UINT8, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!maskarr)
      goto error;

    if(PyArray_Size(maskarr) < PyArray_Size(maparr)) {
      PyErr_SetString(PyExc_IndexError, "mask too small");
      goto error;
    }
  }

  outarr = PyArray_SimpleNew(ndim,
                             dims,
                             NPY_FLOAT);
  if(!outarr)
    goto error;

  if(backremove(PyArray_DATA(maparr),
                maskarr ? PyArray_DATA(maskarr) : NULL,
                PyArray_DATA(outarr),
                dims[DIM_X], dims[DIM_Y], nbsize)) {
    PyErr_SetString(PyExc_MemoryError, "backremove");
    goto error;
  }

  Py_DECREF(maparr);
  Py_XDECREF(maskarr);

  return(PyArray_Return((PyArrayObject *) outarr));

 error:
  Py_XDECREF(maparr);
  Py_XDECREF(maskarr);
  PyArray_XDECREF_ERR((PyArrayObject *) outarr);

  return(NULL);
}

static PyObject *lfa_skylevel_image (PyObject *self,
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

  PyObject *maparr = NULL;
  PyObject *maskarr = NULL;

#ifdef NAN
  float skylev = NAN;
  float skynoise = NAN;
#else
  float skylev = strtod("NaN", (char *) NULL);
  float skynoise = strtod("NaN", (char *) NULL);
#endif

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

  if(maskarg && maskarg != Py_None) {
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

static PyObject *lfa_ap_phot (PyObject *self,
                              PyObject *args,
                              PyObject *kwds) {
  static char *kwlist[] = { "map",
                            "sky",
                            "xcent",
                            "ycent",
                            "r",
                            "mask",
                            NULL };
  PyObject *maparg;
  PyObject *skyarg;
  PyObject *xcentarg;
  PyObject *ycentarg;
  PyObject *rarg;
  PyObject *maskarg = NULL;

  PyObject *maparr = NULL;
  PyObject *skyarr = NULL;
  PyObject *xcentarr = NULL;
  PyObject *ycentarr = NULL;
  PyObject *rarr = NULL;
  PyObject *maskarr = NULL;

  PyObject *xoutarr = NULL;
  PyObject *youtarr = NULL;
  PyObject *fluxoutarr = NULL;

  int iobj, nobj, nsky, nr;

  float *sky;
  double *xcent;
  double *ycent;
  double *r;
  double *xout;
  double *yout;
  double *fluxout;

  int ndim;
  npy_intp *dims;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OOOOO|O", kwlist,
                                  &maparg,
                                  &skyarg,
                                  &xcentarg,
                                  &ycentarg,
                                  &rarg,
                                  &maskarg))
    return(NULL);

  maparr = PyArray_FROM_OTF(maparg, NPY_FLOAT, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!maparr)
    goto error;

  ndim = PyArray_NDIM(maparr);
  dims = PyArray_DIMS(maparr);

  if(ndim != 2) {
    PyErr_SetString(PyExc_IndexError, "expected 2-D map");
    goto error;
  }

  skyarr = PyArray_FROM_OTF(skyarg, NPY_FLOAT, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!skyarr)
    goto error;

  xcentarr = PyArray_FROM_OTF(xcentarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!xcentarr)
    goto error;

  ycentarr = PyArray_FROM_OTF(ycentarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!ycentarr)
    goto error;

  rarr = PyArray_FROM_OTF(rarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!rarr)
    goto error;

  nobj = PyArray_Size(xcentarr);

  if(PyArray_Size(ycentarr) != nobj) {
    PyErr_SetString(PyExc_IndexError, "xcent and ycent sizes differ");
    goto error;
  }

  nsky = PyArray_Size(skyarr);

  if(nsky != 1 && nsky != nobj) {
    PyErr_SetString(PyExc_IndexError, "number of skies differs from objects");
    goto error;
  }

  nr = PyArray_Size(rarr);

  if(nr != 1 && nr != nobj) {
    PyErr_SetString(PyExc_IndexError, "number of skies differs from objects");
    goto error;
  }

  if(maskarg && maskarg != Py_None) {
    maskarr = PyArray_FROM_OTF(maskarg, NPY_UINT8, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!maskarr)
      goto error;

    if(PyArray_Size(maskarr) < PyArray_Size(maparr)) {
      PyErr_SetString(PyExc_IndexError, "mask too small");
      goto error;
    }
  }

  xoutarr = PyArray_SimpleNew(PyArray_NDIM(xcentarr),
                              PyArray_DIMS(xcentarr),
                              NPY_DOUBLE);
  if(!xoutarr)
    goto error;

  youtarr = PyArray_SimpleNew(PyArray_NDIM(xcentarr),
                              PyArray_DIMS(xcentarr),
                              NPY_DOUBLE);
  if(!youtarr)
    goto error;

  fluxoutarr = PyArray_SimpleNew(PyArray_NDIM(xcentarr),
                                 PyArray_DIMS(xcentarr),
                                 NPY_DOUBLE);
  if(!fluxoutarr)
    goto error;

  sky = PyArray_DATA(skyarr);
  xcent = PyArray_DATA(xcentarr);
  ycent = PyArray_DATA(ycentarr);
  r = PyArray_DATA(rarr);
  xout = PyArray_DATA(xoutarr);
  yout = PyArray_DATA(youtarr);
  fluxout = PyArray_DATA(fluxoutarr);

  for(iobj = 0; iobj < nobj; iobj++)
    ap_phot(PyArray_DATA(maparr),
            maskarr ? PyArray_DATA(maskarr) : NULL,
            dims[DIM_X], dims[DIM_Y],
            nsky == 1 ? sky[0] : sky[iobj],
            xcent[iobj], ycent[iobj],
            nr == 1 ? r[0] : r[iobj],
            &(xout[iobj]), &(yout[iobj]), &(fluxout[iobj]));

  Py_DECREF(maparr);
  Py_DECREF(skyarr);
  Py_DECREF(xcentarr);
  Py_DECREF(ycentarr);
  Py_DECREF(rarr);
  Py_XDECREF(maskarr);

  return(Py_BuildValue("NNN",
                       PyArray_Return((PyArrayObject *) xoutarr),
                       PyArray_Return((PyArrayObject *) youtarr),
                       PyArray_Return((PyArrayObject *) fluxoutarr)));

 error:
  Py_XDECREF(maparr);
  Py_XDECREF(skyarr);
  Py_XDECREF(xcentarr);
  Py_XDECREF(ycentarr);
  Py_XDECREF(rarr);
  Py_XDECREF(maskarr);
  PyArray_XDECREF_ERR((PyArrayObject *) xoutarr);
  PyArray_XDECREF_ERR((PyArrayObject *) youtarr);
  PyArray_XDECREF_ERR((PyArrayObject *) fluxoutarr);

  return(NULL);
}

static PyObject *lfa_ap_skyann (PyObject *self,
                                PyObject *args,
                                PyObject *kwds) {
  static char *kwlist[] = { "map",
                            "xcent",
                            "ycent",
                            "rinn",
                            "rout",
                            "mask",
                            "clip_low",
                            "clip_high",
                            NULL };
  PyObject *maparg;
  PyObject *xcentarg;
  PyObject *ycentarg;
  PyObject *rinnarg;
  PyObject *routarg;
  PyObject *maskarg = NULL;

  PyObject *maparr = NULL;
  PyObject *xcentarr = NULL;
  PyObject *ycentarr = NULL;
  PyObject *rinnarr = NULL;
  PyObject *routarr = NULL;
  PyObject *maskarr = NULL;

  PyObject *skylevarr = NULL;
  PyObject *skyrmsarr = NULL;

  int iobj, nobj, nr;

  double *xcent;
  double *ycent;
  double *rinn;
  double *rout;
  float *skylev;
  float *skyrms;

  int ndim;
  npy_intp *dims;

  float clip_low = -FLT_MAX;
  float clip_high = 3.0;

  int hist[SKYLEVEL_DEFAULT_SIZE], hmin, hmax;
  
  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OOOOO|O", kwlist,
                                  &maparg,
                                  &xcentarg,
                                  &ycentarg,
                                  &rinnarg,
                                  &routarg,
                                  &maskarg,
                                  &clip_low,
                                  &clip_high))
    return(NULL);

  maparr = PyArray_FROM_OTF(maparg, NPY_FLOAT, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!maparr)
    goto error;

  ndim = PyArray_NDIM(maparr);
  dims = PyArray_DIMS(maparr);

  if(ndim != 2) {
    PyErr_SetString(PyExc_IndexError, "expected 2-D map");
    goto error;
  }

  xcentarr = PyArray_FROM_OTF(xcentarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!xcentarr)
    goto error;

  ycentarr = PyArray_FROM_OTF(ycentarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!ycentarr)
    goto error;

  rinnarr = PyArray_FROM_OTF(rinnarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!rinnarr)
    goto error;

  routarr = PyArray_FROM_OTF(routarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!routarr)
    goto error;

  nobj = PyArray_Size(xcentarr);

  if(PyArray_Size(ycentarr) != nobj) {
    PyErr_SetString(PyExc_IndexError, "xcent and ycent sizes differ");
    goto error;
  }

  nr = PyArray_Size(rinnarr);

  if(nr != 1 && nr != nobj) {
    PyErr_SetString(PyExc_IndexError, "number of skies differs from objects");
    goto error;
  }

  if(nr != PyArray_Size(routarr)) {
    PyErr_SetString(PyExc_IndexError, "rinn and rout sizes differ");
    goto error;
  }
  
  if(maskarg && maskarg != Py_None) {
    maskarr = PyArray_FROM_OTF(maskarg, NPY_UINT8, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!maskarr)
      goto error;

    if(PyArray_Size(maskarr) < PyArray_Size(maparr)) {
      PyErr_SetString(PyExc_IndexError, "mask too small");
      goto error;
    }
  }

  skylevarr = PyArray_SimpleNew(PyArray_NDIM(xcentarr),
                                PyArray_DIMS(xcentarr),
                                NPY_FLOAT);
  if(!skylevarr)
    goto error;

  skyrmsarr = PyArray_SimpleNew(PyArray_NDIM(xcentarr),
                                PyArray_DIMS(xcentarr),
                                NPY_FLOAT);
  if(!skyrmsarr)
    goto error;

  xcent = PyArray_DATA(xcentarr);
  ycent = PyArray_DATA(ycentarr);
  rinn = PyArray_DATA(rinnarr);
  rout = PyArray_DATA(routarr);
  skylev = PyArray_DATA(skylevarr);
  skyrms = PyArray_DATA(skyrmsarr);

  hmin = 0;
  hmax = SKYLEVEL_DEFAULT_ULIM;

  for(iobj = 0; iobj < nobj; iobj++)
    ap_skyann(PyArray_DATA(maparr),
              maskarr ? PyArray_DATA(maskarr) : NULL,
              dims[DIM_X], dims[DIM_Y],
              xcent[iobj], ycent[iobj],
              nr == 1 ? rinn[0] : rinn[iobj],
              nr == 1 ? rout[0] : rout[iobj],
              hist, &hmin, &hmax,
              clip_low, clip_high,
              &(skylev[iobj]), &(skyrms[iobj]));

  Py_DECREF(maparr);
  Py_DECREF(xcentarr);
  Py_DECREF(ycentarr);
  Py_DECREF(rinnarr);
  Py_DECREF(routarr);
  Py_XDECREF(maskarr);

  return(Py_BuildValue("NN",
                       PyArray_Return((PyArrayObject *) skylevarr),
                       PyArray_Return((PyArrayObject *) skyrmsarr)));

 error:
  Py_XDECREF(maparr);
  Py_XDECREF(xcentarr);
  Py_XDECREF(ycentarr);
  Py_XDECREF(rinnarr);
  Py_XDECREF(routarr);
  Py_XDECREF(maskarr);
  PyArray_XDECREF_ERR((PyArrayObject *) skylevarr);
  PyArray_XDECREF_ERR((PyArrayObject *) skyrmsarr);

  return(NULL);
}

static PyArray_Descr *lfa_specfind_build_descr (void) {
  PyObject *names = (PyObject *) NULL;
  PyObject *formats = (PyObject *) NULL;
  PyObject *offsets = (PyObject *) NULL;
  PyObject *d = (PyObject *) NULL;
  PyArray_Descr *dtype = (PyArray_Descr *) NULL;
  PyObject *o;

  struct {
    char *name;
    char *format;
    size_t offset;
  } members[] = {
    { "x", "f8", offsetof(struct specfind_line, x) },
    { "ziso", "f8", offsetof(struct specfind_line, ziso) },
    { "zpk", "f8", offsetof(struct specfind_line, zpk) },
    { "xl", "i4", offsetof(struct specfind_line, xl) },
    { "xh", "i4", offsetof(struct specfind_line, xh) },
    { "flags", "B", offsetof(struct specfind_line, flags) },
  };

  int imemb, nmemb;

  /* Create lists */
  nmemb = sizeof(members) / sizeof(members[0]);

  names = PyList_New(nmemb);
  formats = PyList_New(nmemb);
  offsets = PyList_New(nmemb);
  if(!names || !formats || !offsets)
    goto error;

  /* Populate */
  for(imemb = 0; imemb < nmemb; imemb++) {
#if PY_MAJOR_VERSION >= 3
    o = PyUnicode_FromString(members[imemb].name);
#else
    o = PyString_FromString(members[imemb].name);
#endif
    if(!o)
      goto error;

    PyList_SetItem(names, imemb, o);  /* list takes over our ref to o */

#if PY_MAJOR_VERSION >= 3
    o = PyUnicode_FromString(members[imemb].format);
#else
    o = PyString_FromString(members[imemb].format);
#endif
    if(!o)
      goto error;

    PyList_SetItem(formats, imemb, o);

#if PY_MAJOR_VERSION >= 3
    o = PyLong_FromSize_t(members[imemb].offset);
#else
    o = PyInt_FromSize_t(members[imemb].offset);
#endif
    if(!o)
      goto error;

    PyList_SetItem(offsets, imemb, o);
  }

  /* Create dictionary */
  d = PyDict_New();
  if(!d)
    goto error;

  PyDict_SetItemString(d, "names", names);  /* does not take over our ref */
  Py_DECREF(names);                         /* so we have to do this */

  PyDict_SetItemString(d, "formats", formats);
  Py_DECREF(formats);

  PyDict_SetItemString(d, "offsets", offsets);
  Py_DECREF(offsets);

#if PY_MAJOR_VERSION >= 3
  o = PyLong_FromSize_t(sizeof(struct specfind_line));
#else
  o = PyInt_FromSize_t(sizeof(struct specfind_line));
#endif
  if(!o)
    goto error;

  PyDict_SetItemString(d, "itemsize", o);
  Py_DECREF(o);

  /* Convert to descriptor */
  if(PyArray_DescrConverter(d, &dtype) != NPY_SUCCEED)
    goto error;

  Py_DECREF(d);

  return(dtype);

 error:
  Py_XDECREF(names);
  Py_XDECREF(formats);
  Py_XDECREF(offsets);
  Py_XDECREF(d);

  return(NULL);
}

/* In order to avoid leaking source lists when we wrap them in numpy
   arrays, we need to create an object to "own" the source list. */

struct lfa_specfind_list_object {
  PyObject_HEAD
  struct specfind_line *list;
  int nlist;
};

static PyObject *lfa_specfind_list_new (PyTypeObject *type,
                                         PyObject *args,
                                         PyObject *kwds) {
  struct lfa_specfind_list_object *self = NULL;

  /* Allocate */
  self = (struct lfa_specfind_list_object *) type->tp_alloc(type, 0);
  if(self) {
    self->list = NULL;
    self->nlist = 0;
  }

  return((PyObject *) self);
}

static void lfa_specfind_list_dealloc (struct lfa_specfind_list_object *self) {
  if(self->list) {
    free(self->list);
    self->list = NULL;
    self->nlist = 0;
  }

  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyTypeObject lfa_specfind_list_type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "lfa.ap.source_list",                      /* tp_name */
  sizeof(struct lfa_specfind_list_object),  /* tp_basicsize */
  0,                                         /* tp_itemsize */
  (destructor) lfa_specfind_list_dealloc,   /* tp_dealloc */
  0,                                         /* tp_print */
  0,                                         /* tp_getattr */
  0,                                         /* tp_setattr */
  0,                                         /* tp_compare */
  0,                                         /* tp_repr */
  0,                                         /* tp_as_number */
  0,                                         /* tp_as_sequence */
  0,                                         /* tp_as_mapping */
  0,                                         /* tp_hash */
  0,                                         /* tp_call */
  0,                                         /* tp_str */
  0,                                         /* tp_getattro */
  0,                                         /* tp_setattro */
  0,                                         /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,  /* tp_flags */
  "ap source list",                          /* tp_doc */
  0,                                         /* tp_traverse */
  0,                                         /* tp_clear */
  0,                                         /* tp_richcompare */
  0,                                         /* tp_weaklistoffset */
  0,                                         /* tp_iter */
  0,                                         /* tp_iternext */
  0,                                         /* tp_methods */
  0,                                         /* tp_members */
  0,                                         /* tp_getset */
  0,                                         /* tp_base */
  0,                                         /* tp_dict */
  0,                                         /* tp_descr_get */
  0,                                         /* tp_descr_set */
  0,                                         /* tp_dictoffset */
  0,                                         /* tp_init */
  0,                                         /* tp_alloc */
  lfa_specfind_list_new                     /* tp_new */
};

static PyObject *lfa_specfind (PyObject *self,
                               PyObject *args,
                               PyObject *kwds) {
  static char *kwlist[] = { "line",
                            "minpix",
                            "sky",
                            "thresh",
                            "filtline",
                            "mask",
                            "overlp",
                            NULL };
  PyObject *linearg;
  PyObject *filtlinearg = NULL;
  PyObject *maskarg = NULL;
  int nx, minpix;
  float sky, thresh;

  PyObject *linearr = NULL;
  PyObject *filtlinearr = NULL;
  PyObject *maskarr = NULL;

  int ioverlp = 1;
  
  struct lfa_specfind_list_object *slist = NULL;
  npy_intp outdim[1] = { 0 };
  PyObject *result = NULL;

  PyArray_Descr *dtype = NULL;
  
  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "Oiff|OOi", kwlist,
                                  &linearg,
                                  &minpix,
                                  &sky,
                                  &thresh,
                                  &filtlinearg,
                                  &maskarg,
                                  &ioverlp))
    return(NULL);

  linearr = PyArray_FROM_OTF(linearg, NPY_FLOAT, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!linearr)
    goto error;

  nx = PyArray_Size(linearr);
  
  if(filtlinearg && filtlinearg != Py_None) {
    filtlinearr = PyArray_FROM_OTF(filtlinearg, NPY_FLOAT, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!filtlinearr)
      goto error;

    if(PyArray_Size(filtlinearr) < nx) {
      PyErr_SetString(PyExc_IndexError, "filtline too small");
      goto error;
    }
  }

  if(maskarg && maskarg != Py_None) {
    maskarr = PyArray_FROM_OTF(maskarg, NPY_UINT8, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!maskarr)
      goto error;

    if(PyArray_Size(maskarr) < nx) {
      PyErr_SetString(PyExc_IndexError, "mask too small");
      goto error;
    }
  }

  /* Create source list object */
  slist = (struct lfa_specfind_list_object *)
    lfa_specfind_list_type.tp_alloc(&lfa_specfind_list_type, 0);
  if(!slist)
    goto error;

  slist->list = NULL;
  slist->nlist = 0;

  if(specfind(PyArray_DATA(linearr),
              filtlinearr ? PyArray_DATA(filtlinearr) : NULL,
              maskarr ? PyArray_DATA(maskarr) : NULL,
              nx,
              minpix, sky, thresh, ioverlp,
              &(slist->list), &(slist->nlist))) {
    PyErr_SetString(PyExc_RuntimeError, "specfind");
    goto error;
  }

  Py_DECREF(linearr);
  Py_XDECREF(filtlinearr);
  Py_XDECREF(maskarr);

  /* Build descriptor */
  dtype = lfa_specfind_build_descr();
  if(!dtype)
    goto error;
  
  /* Create result array.  This takes over our reference to "dtype".
     The memory buffer for the source list is owned by "slist" */
  outdim[0] = slist->nlist;

  result = PyArray_NewFromDescr(&PyArray_Type,
                                dtype,
                                1, outdim, NULL,
                                slist->list, 0,
                                NULL);
  if(!result)
    goto error;

  /* Tell numpy who owns the memory.  This steals our reference to the
     source list, so we don't need to decref it. */
#if defined(NPY_API_VERSION) && defined(NPY_1_7_API_VERSION) && NPY_API_VERSION >= NPY_1_7_API_VERSION
  PyArray_SetBaseObject((PyArrayObject *) result, (PyObject *) slist);
#else
  PyArray_BASE((PyArrayObject *) result) = (PyObject *) slist;
#endif

  /* Finally we can return the result */
  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) result)));

 error:
  Py_XDECREF(linearr);
  Py_XDECREF(filtlinearr);
  Py_XDECREF(maskarr);
  Py_XDECREF((PyObject *) slist);
  Py_XDECREF(dtype);
  PyArray_XDECREF_ERR((PyArrayObject *) result);

  return(NULL);
}

static PyObject *lfa_lrmatch (PyObject *self,
                              PyObject *args,
                              PyObject *kwds) {
  static char *kwlist[] = { "comx",
                            "comy",
                            "comlogrank",
                            "comerr",
                            "refx",
                            "refy",
                            "reflogrank",
                            "referr",
                            "searchrad",
                            "sorted_y",
                            NULL };
  PyObject *comxarg;
  PyObject *comyarg;
  PyObject *comlogrankarg;
  PyObject *comerrarg;
  PyObject *refxarg;
  PyObject *refyarg;
  PyObject *reflogrankarg;
  PyObject *referrarg;
  double searchrad;
  int sorted_y;

  PyObject *comxarr = NULL;
  PyObject *comyarr = NULL;
  PyObject *comlogrankarr = NULL;
  PyObject *comerrarr = NULL;
  PyObject *refxarr = NULL;
  PyObject *refyarr = NULL;
  PyObject *reflogrankarr = NULL;
  PyObject *referrarr = NULL;

  PyObject *best_ref_for_com_arr = NULL;
  PyObject *best_com_for_ref_arr = NULL;

  int ncom, nref, rv;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OOOOOOOOdi", kwlist,
                                  &comxarg,
                                  &comyarg,
                                  &comlogrankarg,
                                  &comerrarg,
                                  &refxarg,
                                  &refyarg,
                                  &reflogrankarg,
                                  &referrarg,
                                  &searchrad,
                                  &sorted_y))
    return(NULL);

  comxarr = PyArray_FROM_OTF(comxarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!comxarr)
    goto error;

  ncom = PyArray_Size(comxarr);

  comyarr = PyArray_FROM_OTF(comyarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!comyarr)
    goto error;

  if(ncom != PyArray_Size(comyarr)) {
    PyErr_SetString(PyExc_IndexError, "comx and comy differ in length");
    goto error;
  }
  
  if(comlogrankarg && comlogrankarg != Py_None) {
    comlogrankarr = PyArray_FROM_OTF(comlogrankarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!comlogrankarr)
      goto error;

    if(ncom != PyArray_Size(comlogrankarr)) {
      PyErr_SetString(PyExc_IndexError, "comx and comlogrank differ in length");
      goto error;
    }
  }

  if(comerrarg && comerrarg != Py_None) {
    comerrarr = PyArray_FROM_OTF(comerrarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!comerrarr)
      goto error;

    if(ncom != PyArray_Size(comerrarr)) {
      PyErr_SetString(PyExc_IndexError, "comx and comerr differ in length");
      goto error;
    }
  }
  
  refxarr = PyArray_FROM_OTF(refxarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!refxarr)
    goto error;

  nref = PyArray_Size(refxarr);

  refyarr = PyArray_FROM_OTF(refyarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!refyarr)
    goto error;

  if(nref != PyArray_Size(refyarr)) {
    PyErr_SetString(PyExc_IndexError, "refx and refy differ in length");
    goto error;
  }
  
  if(reflogrankarg && reflogrankarg != Py_None) {
    reflogrankarr = PyArray_FROM_OTF(reflogrankarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!reflogrankarr)
      goto error;

    if(nref != PyArray_Size(reflogrankarr)) {
      PyErr_SetString(PyExc_IndexError, "refx and reflogrank differ in length");
      goto error;
    }
  }

  if(referrarg && referrarg != Py_None) {
    referrarr = PyArray_FROM_OTF(referrarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!referrarr)
      goto error;

    if(nref != PyArray_Size(referrarr)) {
      PyErr_SetString(PyExc_IndexError, "refx and referr differ in length");
      goto error;
    }
  }

  best_ref_for_com_arr = PyArray_SimpleNew(PyArray_NDIM(comxarr),
                                           PyArray_DIMS(comxarr),
                                           NPY_INT);
  if(!best_ref_for_com_arr)
    goto error;

  best_com_for_ref_arr = PyArray_SimpleNew(PyArray_NDIM(refxarr),
                                           PyArray_DIMS(refxarr),
                                           NPY_INT);
  if(!best_com_for_ref_arr)
    goto error;

  rv = lrmatch(PyArray_DATA(comxarr),
               PyArray_DATA(comyarr),
               comlogrankarr ? PyArray_DATA(comlogrankarr) : NULL,
               comerrarr ? PyArray_DATA(comerrarr) : NULL,
               ncom,
               PyArray_DATA(refxarr),
               PyArray_DATA(refyarr),
               reflogrankarr ? PyArray_DATA(reflogrankarr) : NULL,
               referrarr ? PyArray_DATA(referrarr) : NULL,
               nref,
               searchrad, sorted_y,
               PyArray_DATA(best_ref_for_com_arr),
               PyArray_DATA(best_com_for_ref_arr));
  if(rv < 0) {
    PyErr_SetString(PyExc_MemoryError, "lrmatch");
    goto error;
  }
    
  Py_DECREF(comxarr);
  Py_DECREF(comyarr);
  Py_XDECREF(comlogrankarr);
  Py_XDECREF(comerrarr);
  Py_DECREF(refxarr);
  Py_DECREF(refyarr);
  Py_XDECREF(reflogrankarr);
  Py_XDECREF(referrarr);

  return(Py_BuildValue("NN",
                       PyArray_Return((PyArrayObject *) best_ref_for_com_arr),
                       PyArray_Return((PyArrayObject *) best_com_for_ref_arr)));

 error:
  Py_XDECREF(comxarr);
  Py_XDECREF(comyarr);
  Py_XDECREF(comlogrankarr);
  Py_XDECREF(comerrarr);
  Py_XDECREF(refxarr);
  Py_XDECREF(refyarr);
  Py_XDECREF(reflogrankarr);
  Py_XDECREF(referrarr);
  PyArray_XDECREF_ERR((PyArrayObject *) best_ref_for_com_arr);
  PyArray_XDECREF_ERR((PyArrayObject *) best_com_for_ref_arr);

  return(NULL);
}

static PyObject *lfa_fsgp_kern_sho (PyObject *self,
                                    PyObject *args,
                                    PyObject *kwds) {
  static char *kwlist[] = { "s0",
                            "w0",
                            "q",
                            NULL };

  double s0, w0, q;
  double ktmp[8];
  int nkern;

  PyObject *kernarr = NULL;
  npy_intp outdim[2] = { 0, 4 };

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "ddd", kwlist,
                                  &s0, &w0, &q))
    return(NULL);

  nkern = fsgp_kern_sho(s0, w0, q, ktmp);

  outdim[0] = nkern;

  kernarr = PyArray_SimpleNew(2, outdim, NPY_DOUBLE);
  if(!kernarr)
    goto error;

  memcpy(PyArray_DATA(kernarr), ktmp, PyArray_Size(kernarr) * sizeof(double));

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) kernarr)));

 error:
  PyArray_XDECREF_ERR(kernarr);

  return(NULL);
}

static PyObject *lfa_fsgp_kern_matern (PyObject *self,
                                       PyObject *args,
                                       PyObject *kwds) {
  static char *kwlist[] = { "var",
                            "rho",
                            "eps",
                            NULL };

  double var, rho, eps = 1.0e-3;

  PyObject *kernarr = NULL;
  npy_intp outdim[2] = { 1, 4 };

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "dd|d", kwlist,
                                  &var, &rho, &eps))
    return(NULL);

  kernarr = PyArray_SimpleNew(2, outdim, NPY_DOUBLE);
  if(!kernarr)
    goto error;

  (void) fsgp_kern_matern(var, rho, eps, PyArray_DATA(kernarr));

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) kernarr)));

 error:
  PyArray_XDECREF_ERR(kernarr);

  return(NULL);
}

static PyObject *lfa_fsgp_new (PyTypeObject *type,
                               PyObject *args,
                               PyObject *kwds) {
  struct lfa_fsgp_object *self = NULL;

  /* Allocate */
  self = (struct lfa_fsgp_object *) type->tp_alloc(type, 0);
  return((PyObject *) self);
}

static void lfa_fsgp_dealloc (struct lfa_fsgp_object *self) {
  fsgp_free(&(self->fac));
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *lfa_fsgp_apply (struct lfa_fsgp_object *self,
                                 PyObject *args,
                                 PyObject *kwds) {
  static char *kwlist[] = { "y",
                            NULL };
  PyObject *yarg;
  PyObject *yarr = NULL;
  int nrhs = 1, ndim;
  npy_intp *dims;

  PyObject *zarr = NULL;

  double *fg = NULL;
  double *y, *z, *thisy, *thisz;
  int irhs;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O", kwlist,
                                  &yarg))
    return(NULL);
    
  yarr = PyArray_FROM_OTF(yarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!yarr)
    goto error;

  ndim = PyArray_NDIM(yarr);
  dims = PyArray_DIMS(yarr);

  if(ndim == 1) {
    nrhs = 1;
    if(dims[0] != self->fac.ndp) {
      PyErr_SetString(PyExc_IndexError,
                      "array 'y' dimension != number of data points");
      goto error;
    }
  }
  else if(ndim == 2) {
    nrhs = dims[0];
    if(dims[1] != self->fac.ndp) {
      PyErr_SetString(PyExc_IndexError,
                      "array 'y' second dimension != number of data points");
      goto error;
    }
  }
  else {
    PyErr_SetString(PyExc_IndexError,
                    "don't understand shape of array 'y'");
    goto error;
  }

  zarr = PyArray_SimpleNew(ndim, dims, NPY_DOUBLE);
  if(!zarr)
    goto error;

  /* Allocate workspace */
  fg = (double *) malloc(self->fac.nckern * sizeof(double));
  if(!fg) {
    PyErr_SetString(PyExc_MemoryError, "malloc");
    goto error;
  }

  y = PyArray_DATA(yarr);
  z = PyArray_DATA(zarr);

  for(irhs = 0; irhs < nrhs; irhs++) {
    thisy = y + irhs * self->fac.ndp;
    thisz = z + irhs * self->fac.ndp;

    fsgp_apply(&(self->fac), fg, thisy, thisz);
  }

  free((void *) fg);
  fg = NULL;

  Py_DECREF(yarr);

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) zarr)));

error:
  Py_XDECREF(yarr);
  PyArray_XDECREF_ERR(zarr);

  if(fg)
    free(fg);

  return(NULL);
}

static int lfa_fsgp_compute (struct lfa_fsgp_object *self,
                             PyObject *args,
                             PyObject *kwds) {
  static char *kwlist[] = { "kern",
                            "t",
                            "yerr",
                            NULL };
  PyObject *kernarg;
  PyObject *kernarr = NULL;
  PyObject *targ;
  PyObject *tarr = NULL;
  PyObject *yerrarg;
  PyObject *yerrarr = NULL;

  int nkern = 1, ndim;
  npy_intp *dims;
  int ndp;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "OOO", kwlist,
                                  &kernarg, &targ, &yerrarg))
    return(-1);
    
  kernarr = PyArray_FROM_OTF(kernarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!kernarr)
    goto error;

  ndim = PyArray_NDIM(kernarr);
  dims = PyArray_DIMS(kernarr);

  if(ndim == 1) {
    nkern = 1;
    if(dims[0] != 4) {
      PyErr_SetString(PyExc_IndexError,
                      "array 'kern' dimension != 4");
      goto error;
    }
  }
  else if(ndim == 2) {
    nkern = dims[0];
    if(dims[1] != 4) {
      PyErr_SetString(PyExc_IndexError,
                      "array 'kern' second dimension != 4");
      goto error;
    }
  }
  else {
    PyErr_SetString(PyExc_IndexError,
                    "don't understand shape of array 'kern'");
    goto error;
  }

  tarr = PyArray_FROM_OTF(targ, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!tarr)
    goto error;

  ndp = PyArray_Size(tarr);

  yerrarr = PyArray_FROM_OTF(yerrarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!yerrarr)
    goto error;

  if(PyArray_Size(yerrarr) != ndp) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'yerr' is not the same size as 't'");
    goto error;
  }

  if(fsgp_compute(&(self->fac),
                  PyArray_DATA(kernarr),
                  nkern,
                  PyArray_DATA(tarr),
                  PyArray_DATA(yerrarr),
                  ndp)) {
    PyErr_SetString(PyExc_MemoryError, "fsgp_compute");
    goto error;
  }

  Py_DECREF(kernarr);
  Py_DECREF(tarr);
  Py_DECREF(yerrarr);

  return(0);

error:
  Py_XDECREF(kernarr);
  Py_XDECREF(tarr);
  Py_XDECREF(yerrarr);

  return(-1);
}

static PyObject *lfa_fsgp_logdet (struct lfa_fsgp_object *self,
                                  PyObject *args,
                                  PyObject *kwds) {
  double logdet;
  logdet = fsgp_logdet(&(self->fac));
  return(Py_BuildValue("d", logdet));
}

static PyObject *lfa_fsgp_loglike (struct lfa_fsgp_object *self,
                                   PyObject *args,
                                   PyObject *kwds) {
  static char *kwlist[] = { "y",
                            NULL };
  PyObject *yarg;
  PyObject *yarr = NULL;

  double loglike;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O", kwlist,
                                  &yarg))
    return(NULL);
    
  yarr = PyArray_FROM_OTF(yarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!yarr)
    goto error;

  if(PyArray_Size(yarr) != self->fac.ndp) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'y' size != number of data points");
    goto error;
  }

  if(fsgp_loglike(&(self->fac),
                  PyArray_DATA(yarr),
                  &loglike)) {
    PyErr_SetString(PyExc_MemoryError, "fsgp_loglike");
    goto error;
  }

  Py_DECREF(yarr);

  return(Py_BuildValue("d", loglike));

error:
  Py_XDECREF(yarr);

  return(NULL);
}

static PyObject *lfa_fsgp_predict (struct lfa_fsgp_object *self,
                                   PyObject *args,
                                   PyObject *kwds) {
  static char *kwlist[] = { "y",
                            "tpred",
                            "return_var",
                            NULL };
  PyObject *yarg;
  PyObject *yarr = NULL;
  int npred, ndim;
  npy_intp *dims;

  PyObject *tpredarg = NULL;
  PyObject *tpredarr = NULL;

  PyObject *ypredarr = NULL;
  PyObject *varpredarr = NULL;

  unsigned char return_var = 0;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O|Ob", kwlist,
                                  &yarg, &tpredarg, &return_var))
    return(NULL);
    
  yarr = PyArray_FROM_OTF(yarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!yarr)
    goto error;

  if(PyArray_Size(yarr) != self->fac.ndp) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'y' size != number of data points");
    goto error;
  }

  if(tpredarg && tpredarg != Py_None) {
    tpredarr = PyArray_FROM_OTF(tpredarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!tpredarr)
      goto error;

    npred = PyArray_Size(tpredarr);
    ndim = PyArray_NDIM(tpredarr);
    dims = PyArray_DIMS(tpredarr);
  }
  else {
    npred = PyArray_Size(yarr);
    ndim = PyArray_NDIM(yarr);
    dims = PyArray_DIMS(yarr);
  }

  ypredarr = PyArray_SimpleNew(ndim, dims, NPY_DOUBLE);
  if(!ypredarr)
    goto error;

  if(return_var) {
    varpredarr = PyArray_SimpleNew(ndim, dims, NPY_DOUBLE);
    if(!varpredarr)
      goto error;
  }

  fsgp_predict(&(self->fac),
               PyArray_DATA(yarr),
               tpredarr ? PyArray_DATA(tpredarr) : NULL,
               PyArray_DATA(ypredarr),
               varpredarr ? PyArray_DATA(varpredarr) : NULL,
               npred);

  Py_DECREF(yarr);
  Py_XDECREF(tpredarr);

  if(return_var)
    return(Py_BuildValue("NN",
                         PyArray_Return((PyArrayObject *) ypredarr),
                         PyArray_Return((PyArrayObject *) varpredarr)));
  else
    return(Py_BuildValue("N",
                         PyArray_Return((PyArrayObject *) ypredarr)));

error:
  Py_XDECREF(yarr);
  Py_XDECREF(tpredarr);
  PyArray_XDECREF_ERR(ypredarr);
  PyArray_XDECREF_ERR(varpredarr);

  return(NULL);
}

static PyObject *lfa_fsgp_residual (struct lfa_fsgp_object *self,
                                    PyObject *args,
                                    PyObject *kwds) {
  static char *kwlist[] = { "y",
                            NULL };
  PyObject *yarg;
  PyObject *yarr = NULL;
  int nrhs = 1, ndim;
  npy_intp *dims;

  PyObject *zarr = NULL;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O", kwlist,
                                  &yarg))
    return(NULL);
    
  yarr = PyArray_FROM_OTF(yarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!yarr)
    goto error;

  ndim = PyArray_NDIM(yarr);
  dims = PyArray_DIMS(yarr);

  if(ndim == 1) {
    nrhs = 1;
    if(dims[0] != self->fac.ndp) {
      PyErr_SetString(PyExc_IndexError,
                      "array 'y' dimension != number of data points");
      goto error;
    }
  }
  else if(ndim == 2) {
    nrhs = dims[0];
    if(dims[1] != self->fac.ndp) {
      PyErr_SetString(PyExc_IndexError,
                      "array 'y' second dimension != number of data points");
      goto error;
    }
  }
  else {
    PyErr_SetString(PyExc_IndexError,
                    "don't understand shape of array 'y'");
    goto error;
  }

  zarr = PyArray_SimpleNew(ndim, dims, NPY_DOUBLE);
  if(!zarr)
    goto error;

  fsgp_residual(&(self->fac),
                PyArray_DATA(yarr),
                PyArray_DATA(zarr),
                nrhs);

  Py_DECREF(yarr);

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) zarr)));

error:
  Py_XDECREF(yarr);
  PyArray_XDECREF_ERR(zarr);

  return(NULL);
}

static PyObject *lfa_fsgp_sample (struct lfa_fsgp_object *self,
                                  PyObject *args,
                                  PyObject *kwds) {
  static char *kwlist[] = { "q",
                            NULL };
  PyObject *qarg;
  PyObject *qarr = NULL;
  int ndim;
  npy_intp *dims;

  PyObject *yarr = NULL;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, kwds,
                                  "O", kwlist,
                                  &qarg))
    return(NULL);
    
  qarr = PyArray_FROM_OTF(qarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!qarr)
    goto error;

  if(PyArray_Size(qarr) != self->fac.ndp) {
    PyErr_SetString(PyExc_IndexError,
                    "array 'q' size != number of data points");
    goto error;
  }

  ndim = PyArray_NDIM(qarr);
  dims = PyArray_DIMS(qarr);

  yarr = PyArray_SimpleNew(ndim, dims, NPY_DOUBLE);
  if(!yarr)
    goto error;

  if(fsgp_sample(&(self->fac),
                 PyArray_DATA(qarr),
                 PyArray_DATA(yarr))) {
    PyErr_SetString(PyExc_MemoryError, "fsgp_sample");
    goto error;
  }

  Py_DECREF(qarr);

  return(Py_BuildValue("N",
                       PyArray_Return((PyArrayObject *) yarr)));

error:
  Py_XDECREF(qarr);
  PyArray_XDECREF_ERR(yarr);

  return(NULL);
}

static PyMemberDef lfa_fsgp_members[] = {
  { NULL, 0, 0, 0, NULL }
};

static PyMethodDef lfa_fsgp_methods[] = {
  { "apply", (PyCFunction) lfa_fsgp_apply,
    METH_VARARGS | METH_KEYWORDS,
    "z = apply(y)" },
  { "logdet", (PyCFunction) lfa_fsgp_logdet,
    METH_VARARGS | METH_KEYWORDS,
    "logdet = logdet()" },
  { "loglike", (PyCFunction) lfa_fsgp_loglike,
    METH_VARARGS | METH_KEYWORDS,
    "loglike = loglike(y)" },
  { "predict", (PyCFunction) lfa_fsgp_predict,
    METH_VARARGS | METH_KEYWORDS,
    "z = predict(y, tpred=None, return_var=False)" },
  { "residual", (PyCFunction) lfa_fsgp_residual,
    METH_VARARGS | METH_KEYWORDS,
    "z = residual(y)" },
  { "sample", (PyCFunction) lfa_fsgp_sample,
    METH_VARARGS | METH_KEYWORDS,
    "y = sample(q)" },
  { NULL, NULL, 0, NULL }
};

static PyTypeObject lfa_fsgp_type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "lfa.fsgp",                          /* tp_name */
  sizeof(struct lfa_fsgp_object),      /* tp_basicsize */
  0,                                   /* tp_itemsize */
  (destructor) lfa_fsgp_dealloc,       /* tp_dealloc */
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
  "fsgp structure",                    /* tp_doc */
  0,                                   /* tp_traverse */
  0,                                   /* tp_clear */
  0,                                   /* tp_richcompare */
  0,                                   /* tp_weaklistoffset */
  0,                                   /* tp_iter */
  0,                                   /* tp_iternext */
  lfa_fsgp_methods,                    /* tp_methods */
  lfa_fsgp_members,                    /* tp_members */
  0,                                   /* tp_getset */
  0,                                   /* tp_base */
  0,                                   /* tp_dict */
  0,                                   /* tp_descr_get */
  0,                                   /* tp_descr_set */
  0,                                   /* tp_dictoffset */
  (initproc) lfa_fsgp_compute,         /* tp_init */
  0,                                   /* tp_alloc */
  lfa_fsgp_new                         /* tp_new */
};

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
  { "ranorm", (PyCFunction) lfa_ranorm,
    METH_VARARGS | METH_KEYWORDS,
    "angle = ranorm(angle)" },
  { "range", (PyCFunction) lfa_range,
    METH_VARARGS | METH_KEYWORDS,
    "angle = range(angle)" },  
  { "date2mjd", (PyCFunction) lfa_date2mjd,
    METH_VARARGS | METH_KEYWORDS,
    "mjd = date2mjd(yr, mn, dy)" },
  { "mjd2date", (PyCFunction) lfa_mjd2date,
    METH_VARARGS | METH_KEYWORDS,
    "yr, mn, dy = mjd2date(mjd)" },
  { "geoc", (PyCFunction) lfa_geoc,
    METH_VARARGS | METH_KEYWORDS,
    "u, z = geoc(sinphi, cosphi, height)" },
  { "mount_ab2rp", (PyCFunction) lfa_mount_ab2rp,
    METH_VARARGS | METH_KEYWORDS,
    "pos, r, p = mount_rp2pos(aim, bore, snp=0, cnp=1, flip=0, daimdt=None)" },
  { "mount_pa", (PyCFunction) lfa_mount_pa,
    METH_VARARGS | METH_KEYWORDS,
    "a = mount_pa(aimp, bore, pos)" },
  { "mount_rp2pos", (PyCFunction) lfa_mount_rp2pos,
    METH_VARARGS | METH_KEYWORDS,
    "pos = mount_rp2pos(r, p, snp=0, cnp=1)" },
  { "pfb_matrix", (PyCFunction) lfa_pfb_matrix,
    METH_VARARGS | METH_KEYWORDS,
    "mat = pfb_matrix(mjd_tt)" },
  { "refract_const", (PyCFunction) lfa_refract_const,
    METH_VARARGS | METH_KEYWORDS,
    "refco = refract_const(temperat, humidity, pressure, wavelength, height)" },
  { "refract_corr", (PyCFunction) lfa_refract_corr,
    METH_VARARGS | METH_KEYWORDS,
    "refr, deriv = refract_corr(refco, tanz)" },
  { "refract_vec", (PyCFunction) lfa_refract_vec,
    METH_VARARGS | METH_KEYWORDS,
    "vo = refract_vec(refco, vi, unref=0)" },
  { "dplate", (PyCFunction) lfa_dplate,
    METH_VARARGS | METH_KEYWORDS,
    "tr = dplate(comx, comy, refx, refy, wt=None, ncoeff=6)" },
  { "kepler", (PyCFunction) lfa_kepler,
    METH_VARARGS | METH_KEYWORDS,
    "ea = kepler(ma, ecc)" },
  { "pixovcirc", (PyCFunction) lfa_pixovcirc,
    METH_VARARGS | METH_KEYWORDS,
    "fract = pixovcirc(x, y, r)" },
  { "backremove", (PyCFunction) lfa_backremove,
    METH_VARARGS | METH_KEYWORDS,
    "mapout = backremove(mapin, mask=None, nbsize=64)" },
  { "skylevel_image", (PyCFunction) lfa_skylevel_image,
    METH_VARARGS | METH_KEYWORDS,
    "skylev, skynoise = skylevel_image(map, mask=None, clip_low=-FLT_MAX, clip_high=3)" },
  { "ap_phot", (PyCFunction) lfa_ap_phot,
    METH_VARARGS | METH_KEYWORDS,
    "flux = ap_phot(map, skylev, xcent, ycent, r, mask=None)" },
  { "ap_skyann", (PyCFunction) lfa_ap_skyann,
    METH_VARARGS | METH_KEYWORDS,
    "skylev, skyrms = ap_skyann(map, xcent, ycent, rinn, rout, mask=None, clip_low=-FLT_MAX, clip_high=3)" },
  { "specfind", (PyCFunction) lfa_specfind,
    METH_VARARGS | METH_KEYWORDS,
    "sources = specfind(line, minpix, sky, thresh, filtline=None, mask=None, overlp=1)" },
  { "lrmatch", (PyCFunction) lfa_lrmatch,
    METH_VARARGS | METH_KEYWORDS,
    "best_ref_for_com, best_com_for_ref = lrmatch(comx, comy, comlogrank, comerr, refx, refy, reflogrank, referr, searchrad, sorted_y)" },
  { "fsgp_kern_sho", (PyCFunction) lfa_fsgp_kern_sho,
    METH_VARARGS | METH_KEYWORDS,
    "kern = fsgp_kern_sho(s0, w0, q)" },
  { "fsgp_kern_matern", (PyCFunction) lfa_fsgp_kern_matern,
    METH_VARARGS | METH_KEYWORDS,
    "kern = fsgp_kern_matern(var, rho, eps=1.0e-3)" },
  { NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef lfa_mod = {
  PyModuleDef_HEAD_INIT,
  "lfa",        /* m_name */
  NULL,         /* m_doc */
  -1,           /* m_size */
  lfa_methods,  /* m_methods */
  NULL,         /* m_reload */
  NULL,         /* m_traverse */
  NULL,         /* m_clear */
  NULL,         /* m_free */
};
#endif

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_lfa (void) {
#else
PyMODINIT_FUNC initlfa (void) {
#endif
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
#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&lfa_mod);
#else
  m = Py_InitModule("lfa", lfa_methods);
#endif
  if(!m)
    goto error;

  /* Import numpy */
  import_array();

  /* Create constants */
  nc = sizeof(iconst) / sizeof(iconst[0]);

  for(c = 0; c < nc; c++) {
#if PY_MAJOR_VERSION >= 3
    o = PyLong_FromLong(iconst[c].value);
#else
    o = PyInt_FromLong(iconst[c].value);
#endif
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
    goto error;

  Py_INCREF(&lfa_source_type);
  PyModule_AddObject(m, "source", (PyObject *) &lfa_source_type);

  if(PyType_Ready(&lfa_observer_type) < 0)
    goto error;

  Py_INCREF(&lfa_observer_type);
  PyModule_AddObject(m, "observer", (PyObject *) &lfa_observer_type);

  if(PyType_Ready(&lfa_ap_source_list_type) < 0)
    goto error;

  if(PyType_Ready(&lfa_ap_type) < 0)
    goto error;

  Py_INCREF(&lfa_ap_type);
  PyModule_AddObject(m, "ap", (PyObject *) &lfa_ap_type);

  if(PyType_Ready(&lfa_specfind_list_type) < 0)
    goto error;

  if(PyType_Ready(&lfa_fsgp_type) < 0)
    goto error;

  Py_INCREF(&lfa_fsgp_type);
  PyModule_AddObject(m, "fsgp", (PyObject *) &lfa_fsgp_type);

#if PY_MAJOR_VERSION >= 3
  return(m);

 error:
  return(NULL);
#else
 error:
  return;
#endif
}
