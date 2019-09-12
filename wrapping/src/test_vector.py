#!/usr/bin/env python

#
import sys
import os

import openmeeg_path_setter as ps

import openmeeg as om
import numpy as np

# vector mapping

# om.Vector -> np.array
V1 = om.Vector(3)
V1.set(2.0)
print("V1 =", V1, "\n")

W1 = V1.array()
print("W1 =", W1)

# np.array -> om.Vector
W2 = np.array( [1.0, 2.0, 3.0])
print("W2 of", W2.__class__)
print("W2 =", W2, "\n")

V2 = om.Vector(W2, om.DEEP_COPY)
V2.info()
