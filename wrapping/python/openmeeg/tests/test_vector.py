import numpy as np
import pytest

import openmeeg as om
import openmeeg._openmeeg_wrapper as _omc


def test_vector():
    # _omc.Vector -> np.array
    V1 = _omc.Vector(3)
    V1.set(2.0)
    print("V1 =", V1, V1.info(), "\n")

    W1 = V1.array()
    print("W1 =", W1)
    np.testing.assert_equal(W1, 2 * np.ones(3))

    # np.array -> _omc.Vector
    W2 = np.array([1.0, 2.0, 3.0])
    print("W2 of", W2.__class__)
    print("W2 =", W2, "\n")

    V1p5 = _omc.Vector(_omc.Vector(W2), _omc.DEEP_COPY)
    V1p5.info()
    V2 = _omc.Vector(W2, _omc.DEEP_COPY)  # typemap
    print("V2 of", V2.__class__)
    V2.info()

    V3 = _omc.Vector(V2)
    print("V3 of", V3.__class__)
    V3.info()

    M = om.Matrix(2, 3)
    M.setlin(0, V2)
    M.set(3)
    M.info()
    print("M of", M.__class__)

    print("V3 of", V3.__class__)
    V3 = _omc.Vector(M)
    print("V3 of", V3.__class__)
    V3.info()

    for i in range(3):
        assert W1[i] == V1.value(i)
    print("conversion between OpenMEEG:Vector <> numpy.array is OK")

    # degenerate cases
    with pytest.raises(TypeError, match="Input object is neither.*Vector"):
        _omc.Vector("foo")
    vec = _omc.Vector(3)
    with pytest.raises(IndexError, match="Index out of range"):
        vec.value(3)

    # A 2D array must raise a nice exception, not crash (gh-584): this goes
    # through the same Vector& typemap used e.g. by om.Sensors' weights and
    # radii arguments.
    with pytest.raises(ValueError, match="1 dimensional"):
        _omc.Vector(np.zeros((2, 2)), _omc.DEEP_COPY)


def test_vector_non_contiguous():
    """Read non-contiguous input arrays correctly, not silently misread (gh-584)."""
    padded = np.array([[1.0, -1.0], [2.0, -1.0], [3.0, -1.0]])
    strided = padded[:, 0]
    assert not strided.flags["C_CONTIGUOUS"]
    V = _omc.Vector(strided, _omc.DEEP_COPY)
    np.testing.assert_array_equal(V.array(), strided)
