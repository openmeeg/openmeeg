#!/usr/bin/env python

#
import sys
import os
import numpy

# import locally built openmeeg
topdir = os.getcwd()
omdir  = os.path.join(topdir, "build/wrapping/src/")
if not os.path.exists( os.path.join(omdir , "_openmeeg.so")):
    print("Go to openmeeg topdir before launching this script")
    exit(-1)

sys.path.append(omdir)
import openmeeg as om

#
fileskel = os.path.join(topdir, "data/Head1/Head1.")

# dipole
dipoles = om.Matrix()
dipoles.load( fileskel + "dip")
D = om.asarray(dipoles)
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
