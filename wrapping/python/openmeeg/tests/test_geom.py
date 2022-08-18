#!/usr/bin/env python

import numpy as np
import pytest
import openmeeg as om


def test_geometry():
    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    triangles = np.array([[1, 2, 3], [2, 3, 0]])

    g = om.Geometry()
    mesh = om.Mesh(vertices, triangles, "test", g)

    assert mesh.geometry().check(mesh)

    with pytest.raises(TypeError, match="should be an array"):
        g.add_vertices("foo")
    with pytest.raises(ValueError, match="Vertices.*cannot be converted"):
        g.add_vertices(np.array([1j]))
    with pytest.raises(ValueError, match="Vertices.*2 dim.*3 columns"):
        g.add_vertices(np.array([0.0]))
    with pytest.raises(ValueError, match="Vertices.*2 dim.*3 columns"):
        g.add_vertices(np.array([[0.0]]))
    # TODO should be IOError and have a better error message
    with pytest.raises(IOError, match="Unknown foo suffix"):
        om.Geometry("a.foo")
    with pytest.raises(TypeError, match="Argument.*must be a list"):
        om.Geometry(())
    with pytest.raises(TypeError, match="must be a list of lists"):
        om.Geometry([()])
    with pytest.raises(TypeError, match="first entry a non-empty"):
        om.Geometry([[0, np.zeros((1, 3)), 0]])
