import os.path as op
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


def _assert_geometry(g1, g2):
    assert g1.is_nested()
    assert g2.is_nested()
    assert g1.vertices().size() == g2.vertices().size()
    assert g1.domains().size() == g2.domains().size() == 4
    assert g1.__class__ == g2.__class__

    for d1, d2 in zip(g1.domains(), g2.domains()):
        assert d1.conductivity() == d2.conductivity()
        assert d1.boundaries().size() == d2.boundaries().size()
        for b1, b2 in zip(d1.boundaries(), d2.boundaries()):
            assert b1.inside() == b2.inside()
            assert b1.interface().nb_vertices() == b2.interface().nb_vertices()
            assert b1.interface().nb_triangles() == b2.interface().nb_triangles()
            assert b1.__class__ == b2.__class__


def test_make_geometry_head(data_path):
    # Load mesh data to mimic Head1.geom + Head1.cond
    subject_id = "1"
    subject = f"Head{subject_id}"
    dirpath = op.join(data_path, subject)

    # Make sure we handle bad paths gracefully
    with pytest.raises(IOError, match="Cannot open file"):
        om.Mesh(op.join(dirpath, "fake.1.tri"))

    meshes = list()
    for key in ("cortex", "skull", "scalp"):
        meshes.append(om.Mesh(op.join(dirpath, f"{key}.{subject_id}.tri")))

    g1 = om.make_nested_geometry(meshes, conductivity=(1, 0.0125, 1))
    g2 = om.Geometry(
        op.join(dirpath, subject + ".geom"), op.join(dirpath, subject + ".cond")
    )

    _assert_geometry(g1, g2)
