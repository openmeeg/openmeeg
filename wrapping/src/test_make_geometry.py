#!/usr/bin/env python3
#
# WIP: copy this file in wrapping/src directory of the build tree
# then execute it (python3 ./wrapping/src/test_make_geometry
# path_to_data_directory)

import os.path as op
import numpy as np
from optparse import OptionParser

import openmeeg as om


def python_mesh(name, path):
    mesh = om.Mesh(path)
    mesh_vertices = mesh.geometry().vertices()
    vertices = np.array([vertex.array() for vertex in mesh_vertices])
    mesh_triangles = mesh.triangles()
    triangles = np.array([mesh.triangle(triangle).array() for triangle in mesh_triangles])
    return vertices, triangles

data_path = op.dirname(op.abspath(__file__))
parser = OptionParser()
parser.add_option("-p", "--path", dest="data_path",
                  help="path to data folder", metavar="FILE",
                  default=data_path)

options, args = parser.parse_args()
data_path = options.data_path

# Load mesh data to mimic Head1.geom + Head1.cond

subject = "Head1"
dirpath = op.join(data_path, subject)

meshes = dict()
meshes["cortex"] = python_mesh("cortex", op.join(dirpath, "cortex.1.tri"))
meshes["skull"] = python_mesh("skull", op.join(dirpath, "skull.1.tri"))
meshes["scalp"] = python_mesh("scalp", op.join(dirpath, "scalp.1.tri"))

# It should be possible to have multiple oriented meshes per interface. e.g.
# interface1 = [(m1,om.OrientedMesh.Normal), (m2,om.OrientedMesh.Opposite), (m3,om.OrientedMesh.Normal)]
# It should also be possible to have a name added at the beginning of the
# tuple.

interfaces = {
    "interface1": [('cortex', om.OrientedMesh.Normal)],
    "interface2": [('skull', om.OrientedMesh.Normal)],
    "interface3": [('scalp', om.OrientedMesh.Normal)]
}

domains = {
    "Scalp": ([('interface2', om.SimpleDomain.Outside), ('interface3', om.SimpleDomain.Inside)], 1.0),
    "Brain": ([('interface1', om.SimpleDomain.Inside)], 1.0),
    "Air": ([('interface3', om.SimpleDomain.Outside)], 0.0),
    "Skull": ([('interface2', om.SimpleDomain.Inside), ('interface1', om.SimpleDomain.Outside)], 0.0125)
}

g1 = om.make_geometry(meshes, interfaces, domains)
g2 = om.Geometry(op.join(dirpath, subject + ".geom"), op.join(dirpath, subject + ".cond"))

assert g1.is_nested()
assert g2.is_nested()
assert g1.__class__ == g2.__class__
