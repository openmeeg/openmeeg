#!/usr/bin/env python

#
import sys
import os
import numpy

import openmeeg as om

#
fileskel = os.path.join(ops.omdir, "data/Head1/Head1.")

# dipole
dipoles = om.Matrix()
dipoles.load( fileskel + "dip")
D = dipoles.array()
print("D is a" , D.__class__)
print(D)
# Examples of basic linear algebra
print("Determinant of D is equal to: ", numpy.linalg.det(D))


# TODO: sensors.read() using numpy arrays
# TODO: sensors == [ double ]
sensors = om.Sensors()
sensors.load( fileskel + "squids")
# TODO: D = asarray(sensors).copy...
# TODO: sensors_1 = om.Matrix(D)

# TODO: mesh.read() using numpy arrays
# TODO: mesh == [ double ] , [ int ]
#mesh = om.Mesh()
#mesh.load( fileskel + "tri")
# TODO: V = [...]
# TODO: I = [...]
# TODO: mesh_1 = om.Mesh(V, I)

# TODO: geom.read() using raw python
# TODO: conductivity == { string => double}
# TODO: geometry     == { [ mesh ] , { string => [ int ] } }
#geom = om.Geometry()
#geom.read( fileskel + "geom" , fileskel + "cond" )
