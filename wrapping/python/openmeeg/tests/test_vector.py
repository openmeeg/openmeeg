import numpy as np
import pytest
import openmeeg as om


def test_vector():
    # om.Vector -> np.array
    V1 = om.Vector(3)
    V1.set(2.0)
    print("V1 =", V1, V1.info(), "\n")

    W1 = V1.array()
    print("W1 =", W1)
    np.testing.assert_equal(W1, 2 * np.ones(3))

    # np.array -> om.Vector
    W2 = np.array([1.0, 2.0, 3.0])
    print("W2 of", W2.__class__)
    print("W2 =", W2, "\n")

    V2 = om.Vector(W2, om.DEEP_COPY)
    print("V2 of", V2.__class__)
    V2.info()

    V3 = om.Vector(V2)
    print("V3 of", V3.__class__)
    V3.info()

    M = om.Matrix(2, 3)
    M.setlin(0, V2)
    M.set(3)
    M.info()
    print("M of", M.__class__)

    print("V3 of", V3.__class__)
    V3 = om.Vector(M)
    print("V3 of", V3.__class__)
    V3.info()

    for i in range(3):
        assert W1[i] == V1.value(i)
    print("conversion between OpenMEEG:Vector <> numpy.array is OK")

    # degenerate cases
    with pytest.raises(TypeError, match="Input object is neither.*Vector"):
        om.Vector("foo")
    vec = om.Vector(3)
    with pytest.raises(IndexError, match="Index out of range"):
        vec.value(3)
