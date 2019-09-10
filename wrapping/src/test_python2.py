#!/usr/bin/env python

#
import sys
import os

topdir = os.getcwd()
if not topdir.endswith("/openmeeg"):
    print("Go to openmeeg topdir before launching this script")
    exit(-1)

sys.path.append(topdir + "/build/wrapping/src")

import openmeeg as om


dipoles = om.Matrix()
dipoles.load(topdir+"/data/Head1/Head1.dip")
D = om.asarray(dipoles)
