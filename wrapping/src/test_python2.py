#!/usr/bin/env python

#
import sys
import os
topdir = os.getcwd()
print(topdir)
omdir  = os.path.join(topdir, "build/wrapping/src/")
print(omdir)
print(os.path.join(omdir , "_openmeeg.so"))
if not os.path.exists( os.path.join(omdir , "_openmeeg.so")):
    print("Go to openmeeg topdir before launching this script")
    exit(-1)

sys.path.append(omdir)

import openmeeg as om


dipoles = om.Matrix()
dipoles.load(os.path.join(topdir, "data/Head1/Head1.dip"))
D = om.asarray(dipoles)
