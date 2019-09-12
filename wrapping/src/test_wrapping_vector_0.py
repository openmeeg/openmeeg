#!/usr/bin/env python

#
import sys
import os

# import locally built openmeeg
topdir = os.getcwd()
omdir  = os.path.join(topdir, "build/wrapping/src/")
if not os.path.exists( os.path.join(omdir , "_openmeeg.so")):
    print("Go to openmeeg topdir before launching this script")
    exit(-1)

sys.path.append(omdir)

#
import openmeeg as om
import numpy as np

V = om.Vector(3)
V.set(2.0)
print("V of", V.__class__)
print("V =", V, "\n")

W = np.array( [1.0, 2.0, 3.0])
print("W of", W.__class__)
print("W =", W, "\n")

Labels = [ "toto" ]
Positions = np.array([[ 0 , 1 , 2 ]])
Orientations = np.array([[ -1 , -1 , -2 ]])
Weights = np.array([ 0.5 ])
Radii = np.array([ 1 ])

print("s1=om.Sensors(Labels, Positions, Orientations, Weights, Radii)");
r=om.Vector(Radii,om.DEEP_COPY)
r.info()
