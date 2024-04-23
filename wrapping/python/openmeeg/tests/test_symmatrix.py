import numpy as np
from numpy.testing import assert_array_equal

import openmeeg as om


def _assert_equal(matrix, numpy_array):
    assert matrix.nlin() == numpy_array.shape[0]
    assert matrix.ncol() == numpy_array.shape[1]
    for i in range(matrix.nlin()):
        for j in range(matrix.ncol()):
            assert matrix.value(i, j) == numpy_array[i, j]


def test_symmatrix():
    X = om.SymMatrix(2)
    X.setvalue(0, 0, 1.0)
    X.setvalue(0, 1, 2.0)
    X.setvalue(1, 1, 3.0)
    Z = X.array_flat()
    assert isinstance(Z, np.ndarray)
    assert_array_equal(Z, [1.0, 2.0, 3.0])

    n_rows = 3
    X = om.SymMatrix(n_rows)
    Z = X.array_flat()
    assert Z.shape == ((n_rows * (n_rows + 1) // 2),)

    XX = om.SymMatrix(Z)
    assert_array_equal(XX.array_flat(), Z)
