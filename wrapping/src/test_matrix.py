#!/usr/bin/env python
import openmeeg as om
import numpy as np
import os

W = np.array([[0.0, 0.0, 3.0], [0.0, 2.0, 0.0], [1.0, 0.0, 0.0]])
print("W of", W.__class__)
print("W =", W, "\n")

X = om.Matrix(W)
print("X of", X.__class__)
print("X =", X, "\n")

Z = X.array()
print("Z of", Z.__class__)
print("Z =", Z, "\n")

a = np.array([[1, 2, 3], [4, 5, 6]])
print("a=", a)

b = om.Matrix(a)
assert b.nlin() == 2
assert b.ncol() == 3

c = om.Matrix(b)
assert b.nlin() == 2
assert b.ncol() == 3

import random

nlines = random.randrange(10, 20)
ncols = nlines + 2  # test on not squares matric
mat_numpy = np.ndarray(shape=(nlines, ncols), dtype=float, order="F")
for l in range(nlines):
    for c in range(ncols):
        mat_numpy[l][c] = random.randrange(-100, 100)

mat_om = om.Matrix(mat_numpy)

print("dimensions of mat_numpy: ", mat_numpy.shape)

mat_om.info()

# mimic info()

print("First Values of numpy array")
for l in range(5):
    for c in range(5):
        print(mat_numpy[l][c], end=" ")
    print()


error = 0
for l in range(5):
    for c in range(5):
        if mat_numpy[l][c] != mat_om.value(l, c):
            print("matrices differ at:", l, c)
            error = 1

if error == 0:
    print("conversion between OpenMEEG:Matrix <> numpy.ndarray is OK")
