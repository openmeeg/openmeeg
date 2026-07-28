import os

import numpy as np
import pytest

import openmeeg._openmeeg_wrapper as _omc


def check_mesh(name, vertices, triangles, expected_result):
    try:
        mesh = _omc.Mesh(vertices, triangles)
        mesh.info()
    except Exception:
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


def test_mesh_full(data_path):
    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    triangles = np.array([[1, 2, 3], [2, 3, 0]])

    check_mesh("1", vertices, triangles, True)
    check_mesh("2", np.array([0.0, 0.0, 0.0]), triangles, False)
    check_mesh("3", vertices, np.array([0, 1, 2]), False)

    bad_vertices = np.array([0.0, 0.0])
    check_mesh("4", bad_vertices, triangles, False)

    bad_triangles = np.array([1, 2, 3])
    check_mesh("5", vertices, bad_triangles, False)

    triangles = np.array([[1, 2, 3], [2, 3, 0]], dtype=np.uint64)
    check_mesh("6", vertices, triangles, True)

    triangles = np.array([[1, 2, 3], [2, 3, 0]], dtype=np.int64)
    check_mesh("7", vertices, triangles, True)

    triangles = np.array([[1, 2, 3], [2, 3, 0]], dtype=np.int32)
    check_mesh("8", vertices, triangles, True)

    triangles = np.array([[1, 2, 3], [2, 3, 0]], dtype=np.uint32)
    check_mesh("9", vertices, triangles, True)

    # test X -> should be OK

    data_file = os.path.join(data_path, "Head1", "Head1.tri")
    mesh_X = _omc.Mesh()
    mesh_X.load(data_file)
    mesh_X.info()

    # test Y -> redo with np.array()
    V_X = mesh_X.vertices()
    T_X = mesh_X.triangles()
    assert V_X is not None and T_X is not None
    # V_X.torray() should be possible
    # mesh_Y = _omc.Mesh(V_X, T_X)  # XXX fails for now
    # mesh_Y.info()

    # TODO
    # assert mesh_X.nb_vertices()  == mesh_Y.nb_vertices()
    # assert mesh_X.nb_triangles() == mesh_Y.nb_triangles()
    # V_Y = mesh_Y.vertices()
    # T_Y = mesh_Y.triangles()
    # assert_allclose(V_X, V_Y)
    # assert_allclose(T_X, T_Y)


def _triangle_vertex_sets(mesh):
    # Mesh.update() may reorient triangles (reversing/rotating vertex order
    # within a triangle) for consistency, so compare vertex *sets* per
    # triangle rather than the exact index order.
    return sorted(
        frozenset(mesh.triangle(t).array().tolist()) for t in mesh.triangles()
    )


@pytest.mark.parametrize(
    "dtype",
    [
        "<i4",
        ">i4",
        "<u4",
        ">u4",
        "<i8",
        ">i8",
        "<u8",
        ">u8",
    ],
)
def test_mesh_triangles_byteorder(dtype):
    """Read triangle indices identically regardless of dtype/byte order (gh-817)."""
    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    triangles_native = np.array([[1, 2, 3], [2, 3, 0]], dtype=np.int64)
    want = _triangle_vertex_sets(_omc.Mesh(vertices, triangles_native))

    triangles = triangles_native.astype(dtype)
    got = _triangle_vertex_sets(_omc.Mesh(vertices, triangles))
    assert got == want


def test_mesh_triangles_non_contiguous():
    """Non-contiguous (strided) triangle arrays must be read correctly."""
    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    triangles_native = np.array([[1, 2, 3], [2, 3, 0]], dtype=np.int64)
    want = _triangle_vertex_sets(_omc.Mesh(vertices, triangles_native))

    # Add a padding column and slice it back off with basic (view-producing)
    # slicing, so the resulting array is not C-contiguous.
    padded = np.array([[1, 2, 3, 99], [2, 3, 0, 99]], dtype=np.int64)
    triangles = padded[:, :3]
    assert not triangles.flags["C_CONTIGUOUS"]
    got = _triangle_vertex_sets(_omc.Mesh(vertices, triangles))
    assert got == want
