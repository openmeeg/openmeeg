#!/usr/bin/env python
import numpy as np

import openmeeg as om

from os import path as op
from optparse import OptionParser


data_path = op.dirname(op.abspath(__file__))
parser = OptionParser()
parser.add_option(
    "-p",
    "--path",
    dest="data_path",
    help="path to data folder",
    metavar="FILE",
    default=data_path,
)
options, args = parser.parse_args()
data_path = options.data_path

#
subject = "Head1"
file_name_skeleton = op.join(data_path, subject, subject)

# dipole
dipoles = om.Matrix()
dipoles.load(file_name_skeleton + ".dip")
D = dipoles.array()
print("D is a", D.__class__)
print(D)
# Examples of basic linear algebra
print("Determinant of D is equal to: ", np.linalg.det(D))


# TODO: sensors.read() using numpy arrays
# TODO: sensors == [ double ]
sensors = om.Sensors()
sensors.load(file_name_skeleton + ".squids")
# TODO: D = asarray(sensors).copy...
# TODO: sensors_1 = om.Matrix(D)

# TODO: mesh.read() using numpy arrays
# TODO: mesh == [ double ] , [ int ]
# mesh = om.Mesh()
# mesh.load( file_name_skeleton + ".tri")
# TODO: V = [...]
# TODO: I = [...]
# TODO: mesh_1 = om.Mesh(V, I)

# TODO: geom.read() using raw python
# TODO: conductivity == { string => double}
# TODO: geometry     == { [ mesh ] , { string => [ int ] } }
# geom = om.Geometry()
# geom.read( file_name_skeleton + ".geom" , file_name_skeleton + ".cond" )
