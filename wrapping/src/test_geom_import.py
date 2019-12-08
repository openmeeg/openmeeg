#!/usr/bin/env python

import numpy as np

import openmeeg as om

vertices = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                     [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
triangles = np.array([[1, 2, 3], [2, 3, 0]])

mesh = om.Mesh(vertices, triangles)

g = om.Geometry()

assert g.check(mesh)

g.import_meshes([mesh])

assert not g.check(mesh)
