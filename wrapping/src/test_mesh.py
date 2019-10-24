#!/usr/bin/env python

#
import sys
import os
import openmeeg as om
import numpy as np

# test 1 -> should be OK
V_OK = np.array(
    [
        [ 0.0 ,  0.0,  0.0 ],
        [ 1.0 ,  0.0,  0.0 ],
        [ 1.0 ,  1.0,  0.0 ],
        [ 0.0 ,  1.0,  0.0 ]
    ]
)
T_OK = np.array(
    [
        [ 1, 2, 3],
        [ 2, 3, 0]
    ]
)
print("== Test 1")
mesh_1 = om.Mesh(V_OK, T_OK)
mesh_1.info()

# test 2 -> should send a PyError
V_BAD = np.array(
        [ 0.0 ,  0.0,  0.0 ]
)
try:
    print("== Test 2")
    mesh_2 = om.Mesh(V_BAD, T_OK)
    # should not reach this line
    sys.exit(-1)
except:
    print("PyError trapped => OK")


# test 3 -> should send a PyError
T_BAD = np.array(
        [ 0 ,  1,  2 ]
)
try:
    print("== Test 3")
    mesh_3 = om.Mesh(V_OK, T_BAD)
    # should not reach this line
    sys.exit(-1)
except:
    print("PyError trapped => OK")

# test 4 -> should send a PyError
V_BAD = np.array(
    [
        [ 0.0 ,  0.0       ],
        [ 1.0 ,  0.0,  0.0 ],
        [ 1.0 ,  1.0,  0.0 ],
        [ 0.0 ,  1.0,  0.0 ]
    ]
)
try:
    print("== Test 4")
    mesh_4 = om.Mesh(V_BAD, T_OK)
    # should not reach this line
    sys.exit(-1)
except:
    print("PyError trapped => OK")

# test 5 -> should send a PyError
T_BAD = np.array(
    [
        [ 1 ,  2     ],
        [ 0 ,  1,  2 ]
    ]
)
try:
    print("== Test 5")
    mesh_4 = om.Mesh(V_OK, T_BAD)
    # should not reach this line
    sys.exit(-1)
except:
    print("PyError trapped => OK")


sys.exit()

V = np.array(
    [
        [-0.3575,  0,       0.5784,  0.525768,  0,        -0.850628],
        [ 0.3575,  0,       0.5784, -0.525768,  0,        -0.850628],
        [-0.3575,  0,      -0.5784,  0.525768,  0,         0.850628],
        [ 0.3575,  0,      -0.5784, -0.525768,  0,         0.850628],
        [ 0,       0.5784,  0.3575,  0,        -0.850628, -0.525768],
        [ 0,       0.5784, -0.3575,  0,        -0.850628,  0.525768],
        [ 0,       0.68,    0,       0,        -1,         0],
        [ 0.34,    0.5501, -0.2101, -0.50003,  -0.808998,  0.309019],
        [ 0.34,    0.5501,  0.2101, -0.50003,  -0.808998, -0.309019],
        [ 0.5501,  0.2101,  0.34,   -0.808998, -0.309019, -0.50003],
        [ 0.68,    0,       0,      -1,         0,         0],
        [ 0.5501, -0.2101,  0.34,   -0.808998,  0.309019, -0.50003],
        [ 0.5501,  0.2101, -0.34,   -0.808998, -0.309019,  0.50003],
        [ 0.5501, -0.2101, -0.34,   -0.808998,  0.309019,  0.50003],
        [ 0.2101,  0.34,   -0.5501, -0.309019, -0.50003,   0.808998],
        [-0.2101,  0.34,   -0.5501,  0.309019, -0.50003,   0.808998],
        [ 0,       0,      -0.68,    0,         0,         1],
        [-0.2101, -0.34,   -0.5501,  0.309019,  0.50003,   0.808998],
        [ 0.2101, -0.34,   -0.5501, -0.309019,  0.50003,   0.808998],
        [ 0.34,   -0.5501, -0.2101, -0.50003,   0.808998,  0.309019],
        [ 0,      -0.68,    0,       0,         1,         0],
        [ 0.34,   -0.5501,  0.2101, -0.50003,   0.808998, -0.309019],
        [-0.34,   -0.5501, -0.2101,  0.50003,   0.808998,  0.309019],
        [-0.34,   -0.5501,  0.2101,  0.50003,   0.808998, -0.309019],
        [-0.5501, -0.2101,  0.34,    0.808998,  0.309019, -0.50003],
        [-0.2101, -0.34,    0.5501,  0.309019,  0.50003,  -0.808998],
        [ 0.2101, -0.34,    0.5501, -0.309019,  0.50003,  -0.808998],
        [ 0.68,    0,       0.01,    1,         0,         0],
        [-0.5501, -0.2101, -0.34,    0.808998,  0.309019,  0.50003],
        [-0.5501,  0.2101, -0.34,    0.808998, -0.309019,  0.50003],
    ]
)
I = np.array([[1, 3, 4], [2, 4, 5], [10, 3, 5], [21, 4, 5]])
mesh_1 = om.Mesh(V, I)
mesh_1.info()

mesh_f = om.Mesh()
mesh_f.load("/Users/jls/Development/athena/openmeeg/data/Head1/Head1.tri")

# TODO: mesh.read() using numpy arrays
# TODO: mesh == [ double ] , [ int ]
# mesh = om.Mesh()
# mesh.load( fileskel + "tri")
# TODO: V = [...]
# TODO: I = [...]
# TODO: mesh_1 = om.Mesh(V, I)
