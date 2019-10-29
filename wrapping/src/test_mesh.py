#!/usr/bin/env python

#
import sys
import os
import openmeeg as om
import numpy as np

# test 1 -> should be OK
V1_OK = np.array(
    [
        [ 0.0 ,  0.0,  0.0 ],
        [ 1.0 ,  0.0,  0.0 ],
        [ 1.0 ,  1.0,  0.0 ],
        [ 0.0 ,  1.0,  0.0 ]
    ]
)
T1_OK = np.array(
    [
        [ 1, 2, 3],
        [ 2, 3, 0]
    ],
    dtype = np.uint64
)
trap = False
try:
    print("== Test 1")
    mesh_1 = om.Mesh(V1_OK, T1_OK)
    mesh_1.update()
    mesh_1.info()
except:
    # should not reach this line
    trap = True
assert trap == False

# test 1 bis -> should be OK too
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
    ],
    dtype = np.int64
)
trap = False
try:
    print("== Test 1 bis")
    mesh_1 = om.Mesh(V_OK, T_OK)
    mesh_1.update()
    mesh_1.info()
except:
    # should not reach this line
    trap = True
assert trap == False

# test 2 -> should send a PyError
trap = False
V_BAD = np.array(
        [ 0.0 ,  0.0,  0.0 ]
)
try:
    print("== Test 2")
    mesh_2 = om.Mesh(V_BAD, T_OK)
    # should not reach this line
except:
    print("PyError trapped => OK")
    trap = True
assert trap == True

# test 3 -> should send a PyError
trap = False
T_BAD = np.array(
        [ 0 ,  1,  2 ]
)
try:
    print("== Test 3")
    mesh_3 = om.Mesh(V_OK, T_BAD)
    # should not reach this line
except:
    print("PyError trapped => OK")
    trap = True
assert trap == True

# test 4 -> should send a PyError
trap = False
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
except:
    print("PyError trapped => OK")
    trap = True
assert trap == True

# test 5 -> should send a PyError
trap = False
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
except:
    print("PyError trapped => OK")
    trap = True
assert trap == True

# test 6 -> should be OK
test_dir  = os.path.dirname(os.path.abspath(__file__))
data_file = os.path.join( test_dir , "..", "..", "..", "data", "Head1" , "Head1.tri" )
mesh_6 = om.Mesh()
mesh_6.load(data_file)
mesh_6.update()
mesh_6.info()

# test 7 -> redo with np.array()
V6 = mesh_6.vertices()
#T6 = mesh_6.triangles()
#mesh_7 = om.Mesh(V6, T6)
#mesh_7.info()

# TODO
#
# mesh_6.nb_vertices()  == mesh_7.nb_vertices()
# mesh_6.nb_triangles() == mesh_7.nb_triangles()
# V7 = mesh_6.vertices()
# T7 = mesh_6.triangles()
# V6 == V7
# T6 == T7
