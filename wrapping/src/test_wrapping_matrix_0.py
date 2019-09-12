#!/usr/bin/env python

import path_setter as ps
#
import openmeeg as om
import numpy as np
import os

V = om.Matrix(3,3)
V.setvalue(1,1,2)
V.setvalue(2,2,2)
V.setvalue(3,3,4)
print("V of", V.__class__)
print("V =", V, "\n")

W = np.array( [
    [0.0, 0.0, 3.0],
    [0.0, 2.0, 0.0],
    [1.0, 0.0, 0.0]
    ]
)
print("W of", W.__class__)
print("W =", W, "\n")

W = V
print("W of", W.__class__)
print("W =", W, "\n")

#
fileskel = os.path.join(ps.topdir, "data/Head1/Head1.")

S = om.Sensors()
S.load( fileskel + "")
Pos = S.getPosition(1)

S2 = om.Sensors("toto", Pos)
#print("S2 of", S2.__class__)
