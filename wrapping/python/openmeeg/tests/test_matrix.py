import pytest
import numpy as np
from numpy.testing import assert_array_equal
import openmeeg as om


def _assert_equal(matrix, numpy_array):
    assert matrix.nlin() == numpy_array.shape[0]
    assert matrix.ncol() == numpy_array.shape[1]
    for i in range(matrix.nlin()):
        for j in range(matrix.ncol()):
            assert matrix.value(i, j) == numpy_array[i, j]


def test_matrix():
    rng = np.random.RandomState(0)
    W = np.asfortranarray([[0.0, 0.0, 3.0], [0.0, 2.0, 0.0], [1.0, 0.0, 0.0]])

    X = om.Matrix(W)
    Z = X.array()
    assert isinstance(Z, np.ndarray)
    assert_array_equal(Z, W)

    a = np.asfortranarray([[1, 2, 3], [4, 5, 6]])
    b = om.Matrix(a)
    assert b.nlin() == 2
    assert b.ncol() == 3

    for i in range(b.nlin()):
        for j in range(b.ncol()):
            assert b.value(i, j) == a[i, j]

    c = om.Matrix(b)
    assert c.nlin() == 2
    assert c.ncol() == 3

    _assert_equal(c, a)

    nlines = rng.randint(10, 20)
    ncols = nlines + 2  # test on not squared matric

    mat_numpy = np.asfortranarray(rng.randn(nlines, ncols))
    mat_om = om.Matrix(mat_numpy)

    assert (mat_om.nlin(), mat_om.ncol()) == mat_numpy.shape
    mat_om.info()

    _assert_equal(mat_om, mat_numpy)

    # Testing going back from OpenMEEG to numpy
    mat_om2np = mat_om.array()
    assert_array_equal(mat_numpy, mat_om2np)

    with pytest.raises(TypeError, match="can only have 2 dim"):
        om.Matrix(np.zeros((1, 1, 1)))
    with pytest.raises(TypeError, match="Fortran order"):
        om.Matrix(np.zeros((2, 2)))
    with pytest.raises(TypeError, match="Input object must.*OpenMEEG Matrix"):
        om.Matrix([])
    mat = om.Matrix(np.zeros((2, 2)).T)  # okay
    with pytest.raises(IndexError, match="out of range"):
        mat.value(2, 2)
