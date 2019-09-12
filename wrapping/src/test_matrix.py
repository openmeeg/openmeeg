#!/usr/bin/env python

import openmeeg_path_setter as ps

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

X = om.Matrix(W)
print("X of", X.__class__)
print("X =", X, "\n")

Z = X.array()
print("Z of", Z.__class__)
print("Z =", Z, "\n")

a = np.array([[1,2,3],[4,5,6]])
print("a=", a)

b = om.Matrix(a)
assert(b.nlin() == 2)
assert(b.ncol() == 3)

import random
nlines = random.randrange(10,50)
ncols = random.randrange(10,50)
mat_numpy = np.ndarray(shape=(nlines,ncols), dtype=float, order='F')
for l in range(nlines):
    for c in range(ncols):
        print(l, " ", c)
        mat_numpy[l][c] = random.randrange(-100,100)

mat_om = om.Matrix(mat_numpy)

print("dimensions of mat_numpy: ", mat_numpy.shape)
mat_om.info()
