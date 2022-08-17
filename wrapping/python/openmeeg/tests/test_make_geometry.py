#!/usr/bin/env python3

import os.path as op
import numpy as np
import pytest
import openmeeg as om


def test_make_geometry(data_path):
    def python_mesh(path):
        mesh = om.Mesh(path)
        mesh_vertices = mesh.geometry().vertices()
        vertices = np.array([vertex.array() for vertex in mesh_vertices])
        mesh_triangles = mesh.triangles()
        triangles = np.array(
            [mesh.triangle(triangle).array() for triangle in mesh_triangles]
        )
        return vertices, triangles

    # Load mesh data to mimic Head1.geom + Head1.cond
    subject = "Head1"
    dirpath = op.join(data_path, subject)

    # Make sure we handle bad paths gracefully
    with pytest.raises(IOError, match="Cannot open file"):
        om.Mesh(op.join(dirpath, "fake.1.tri"))
    meshes = dict()
    for key in ("cortex", "skull", "scalp"):
        meshes[key] = python_mesh(op.join(dirpath, f"{key}.1.tri"))

    # It should be possible to have multiple oriented meshes per interface. e.g.
    # interface1 = [(m1,om.OrientedMesh.Normal),
    #               (m2,om.OrientedMesh.Opposite),
    #               (m3,om.OrientedMesh.Normal)]
    # It should also be possible to have a name added at the beginning of the
    # tuple.

    interfaces = {
        "interface1": [("cortex", om.OrientedMesh.Normal)],
        "interface2": [("skull", om.OrientedMesh.Normal)],
        "interface3": [("scalp", om.OrientedMesh.Normal)],
    }

    domains = {
        "Scalp": (
            [
                ("interface2", om.SimpleDomain.Outside),
                ("interface3", om.SimpleDomain.Inside),
            ],
            1.0,
        ),
        "Brain": ([("interface1", om.SimpleDomain.Inside)], 1.0),
        "Air": ([("interface3", om.SimpleDomain.Outside)], 0.0),
        "Skull": (
            [
                ("interface2", om.SimpleDomain.Inside),
                ("interface1", om.SimpleDomain.Outside),
            ],
            0.0125,
        ),
    }

    g1 = om.make_geometry(meshes, interfaces, domains)
    g2 = om.Geometry(
        op.join(dirpath, subject + ".geom"), op.join(dirpath, subject + ".cond")
    )

    assert g1.is_nested()
    assert g2.is_nested()
    assert g1.__class__ == g2.__class__
