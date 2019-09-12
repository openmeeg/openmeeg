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

X = om.Vector(W)
print("X of", X.__class__)
print("X =", X, "\n")

Z = X.array()
print("Z of", Z.__class__)
print("Z =", Z, "\n")

Labels = [ "toto" ]
Positions = np.array([ 0 , 1 , 2 ])
Orientations = np.array([ -1 , -1 , -2 ])
Weights = np.array([ 0.5 ])
Radii = np.array([ 1 ])

#S = om.Sensors(Labels, Positions, Orientations, Weights, Radii)
print("om.Sensors(Labels, Positions, Orientations, Weights, Radii)");

#
#fileskel = os.path.join(topdir, "data/Head1/Head1.")

#S = om.Sensors()
#S.load( fileskel + "squids")
#Pos = S.getPosition(1)

#S2 = om.Sensors("toto", Pos)
#print("S2 of", S2.__class__)
