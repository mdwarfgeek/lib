#!/usr/bin/env python

import distutils.core
import distutils.ccompiler
import numpy

np_inc = [ numpy.get_include() ]

# Fetch the list of source files from the Makefile.
# This little parser is rather dumb but should work okay for the way
# the Makefile was written.
srcs = ["pywrap.c"]

with open("Makefile", "r") as fp:
  for line in fp:
    buf = line.strip()
    if buf == "":
      continue

    ll = buf.split("=", 2)
    if len(ll) == 2:
      key = ll[0].strip()
      if key == "SRCS":
        value = ll[1].strip()
        srcs.extend(value.split())

mod = distutils.core.Extension("lfa",
                               include_dirs=np_inc,
                               libraries=["m"],
                               sources=srcs)

distutils.core.setup(name="lfa",
                     version="0.01",
                     description="JMI's library",
                     ext_modules=[mod])

