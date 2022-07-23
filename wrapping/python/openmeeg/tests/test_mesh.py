import os
from os import path as path

import numpy as np
import pytest
import openmeeg as om


@pytest.mark.skipif(os.getenv('OPENMEEG_BAD_PYPY') == '1',
                    reason="bug with PyPy support")
def test_mesh_full(data_path):
    def test_mesh(name, vertices, triangles, expected_result):
        try:
            mesh = om.Mesh(vertices, triangles)
            mesh.info()
        except:
            if expected_result:
                print("Test", name, "--> Failed")
                assert False
            else:
                print("Test", name, "--> Expected failure")
            return
        if not expected_result:
            print("Test", name, "--> Unexpected success")
            assert False
        print("Test", name, "--> Expected success")

    vertices = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                        [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    triangles = np.array([[1, 2, 3], [2, 3, 0]])

    test_mesh("1", vertices, triangles, True)
    test_mesh("2", np.array([0.0, 0.0, 0.0]), triangles, False)
    test_mesh("3", vertices, np.array([0, 1, 2]), False)

    bad_vertices = np.array([0.0, 0.0])
    test_mesh("4", bad_vertices, triangles, False)

    bad_triangles = np.array([1, 2, 3])
    test_mesh("5", vertices, bad_triangles, False)

    triangles = np.array([[1, 2, 3], [2, 3, 0]], dtype=np.uint64)
    test_mesh("6", vertices, triangles, True)

    triangles = np.array([[1, 2, 3], [2, 3, 0]], dtype=np.int64)
    test_mesh("7", vertices, triangles, True)

    # test X -> should be OK

    data_file = path.join(data_path, "Head1", "Head1.tri")
    mesh_X = om.Mesh()
    mesh_X.load(data_file)

    # test Y -> redo with np.array()
    V_Y = mesh_X.vertices()
    # T6 = mesh_6.triangles()
    # mesh_7 = om.Mesh(V6, T6)
    # mesh_7.info()

    # TODO
    #
    # mesh_6.nb_vertices()  == mesh_7.nb_vertices()
    # mesh_6.nb_triangles() == mesh_7.nb_triangles()
    # V7 = mesh_6.vertices()
    # T7 = mesh_6.triangles()
    # V6 == V7
    # T6 == T7
