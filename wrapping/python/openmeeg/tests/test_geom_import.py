#!/usr/bin/env python

import numpy as np

import openmeeg as om


def test_geom_import():
    vertices = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                        [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    triangles = np.array([[1, 2, 3], [2, 3, 0]])

    g = om.Geometry()
    mesh = om.Mesh(vertices, triangles, "test", g)

    assert mesh.geometry().check(mesh)
