#!/usr/bin/env python

#
import sys
import os

omdir  = os.getcwd() + "/build/wrapping/src"

if not os.path.exists( omdir + "_openmeeg.so"):
    print("Go to openmeeg topdir before launching this script")
    exit(-1)

omdir  = topdir + "/build/wrapping/src"
sys.path.append(topdir + "/build/wrapping/src")

import openmeeg as om


dipoles = om.Matrix()
dipoles.load(topdir+"/data/Head1/Head1.dip")
D = om.asarray(dipoles)
D
