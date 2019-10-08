#!/usr/bin/env python
import sys
import os
import numpy as np

import openmeeg as om

# vector mapping

# om.Vector -> np.array
V1 = om.Vector(3)
V1.set(2.0)
print("V1 =", V1, V1.info(), "\n")

W1 = V1.array()
print("W1 =", W1)

# np.array -> om.Vector
W2 = np.array([1.0, 2.0, 3.0])
print("W2 of", W2.__class__)
print("W2 =", W2, "\n")

V2 = om.Vector(W2, om.DEEP_COPY)
print("V2 of", V2.__class__)
V2.info()

V3 = om.Vector(V2)
print("V3 of", V3.__class__)
V3.info()

M = om.Matrix(2, 3)
M.setlin(0, V2)
M.set(3)
M.info()
print("M of", M.__class__)

print("V3 of", V3.__class__)
V3 = om.Vector(M)
print("V3 of", V3.__class__)
V3.info()

#
error = 0
for i in range(3):
    if W1[i] != V1.value(i):
        print("vectors differ at:", i)
        error = 1

if error == 0:
    print("conversion between OpenMEEG:Vector <> numpy.array is OK")
