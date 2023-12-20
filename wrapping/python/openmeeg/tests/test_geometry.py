import os.path as op

import numpy as np
import pytest

import openmeeg as om
from openmeeg._openmeeg_wrapper import Geometry, Mesh


def test_geometry_private():
    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    triangles = np.array([[1, 2, 3], [2, 3, 0]])

    g = Geometry()
    mesh = Mesh(vertices, triangles, "test", g)

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
        Geometry("a.foo")
    with pytest.raises(TypeError, match="Argument.*must be a list"):
        Geometry(())
    with pytest.raises(TypeError, match="must be a list of lists"):
        Geometry([()])
    with pytest.raises(TypeError, match="first entry a non-empty"):
        Geometry([[0, np.zeros((1, 3)), 0]])


def _assert_geometry(g1, g2, n_domains):
    assert g1.is_nested()
    assert g2.is_nested()
    assert g1.vertices().size() == g2.vertices().size()
    assert g1.domains().size() == g2.domains().size() == n_domains
    assert g1.__class__ == g2.__class__

    for d1, d2 in zip(g1.domains(), g2.domains()):
        assert d1.conductivity() == d2.conductivity()
        assert d1.boundaries().size() == d2.boundaries().size()
        for b1, b2 in zip(d1.boundaries(), d2.boundaries()):
            i1, i2 = b1.interface(), b2.interface()
            assert b1.inside() == b2.inside()
            assert i1.nb_vertices() == i2.nb_vertices()
            assert i1.nb_triangles() == i2.nb_triangles()
            assert b1.__class__ == b2.__class__


def test_make_geometry_head(data_path):
    # Load mesh data to mimic Head1.geom + Head1.cond
    subject_id = "1"
    subject = f"Head{subject_id}"
    dirpath = op.join(data_path, subject)

    # Make sure we handle bad paths gracefully
    with pytest.raises(IOError, match="Cannot open file"):
        Mesh(op.join(dirpath, "fake.1.tri"))

    meshes = list()
    for key in ("cortex", "skull", "scalp"):
        meshes.append(Mesh(op.join(dirpath, f"{key}.{subject_id}.tri")))

    # Make a geometry from a 3 layers model
    g1 = om.make_nested_geometry(meshes, conductivity=(1, 0.0125, 1))
    g2 = om.read_geometry(
        op.join(dirpath, subject + ".geom"), op.join(dirpath, subject + ".cond")
    )
    _assert_geometry(g1, g2, n_domains=4)

    # Make a geometry from a 1 layer model
    g1 = om.make_nested_geometry(meshes[:1], conductivity=(1,))
    g2 = om.read_geometry(
        op.join(dirpath, subject + "_1_layer.geom"),
        op.join(dirpath, subject + "_1_layer.cond"),
    )
    _assert_geometry(g1, g2, n_domains=2)
