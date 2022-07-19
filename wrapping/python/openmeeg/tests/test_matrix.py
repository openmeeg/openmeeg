import random
import numpy as np
import openmeeg as om


def test_matrix():
    W = np.asfortranarray([[0.0, 0.0, 3.0], [0.0, 2.0, 0.0], [1.0, 0.0, 0.0]])
    print("W of", W.__class__)
    print("W =", W, "\n")

    X = om.Matrix(W)
    print("X of", X.__class__)
    print("X =", X, "\n")

    Z = X.array()
    print("Z of", Z.__class__)
    print("Z =", Z, "\n")

    a = np.asfortranarray([[1, 2, 3], [4, 5, 6]])
    print("a=", a)

    b = om.Matrix(a)
    assert b.nlin() == 2
    assert b.ncol() == 3

    for i in range(b.nlin()):
        for j in range(b.ncol()):
            assert b.value(i, j) == a[i, j]

    c = om.Matrix(b)
    assert c.nlin() == 2
    assert c.ncol() == 3

    for i in range(c.nlin()):
        for j in range(c.ncol()):
            assert c.value(i, j) == a[i, j]

    nlines = random.randrange(10, 20)
    ncols = nlines + 2  # test on not squared matric

    mat_numpy = np.asfortranarray(np.random.randn(nlines, ncols))
    mat_om = om.Matrix(mat_numpy)

    print("dimensions of mat_numpy: ", mat_numpy.shape)

    assert (mat_om.nlin(), mat_om.ncol()) == mat_numpy.shape

    mat_om.info()

    # mimic info()

    print("First Values of numpy array")
    for l in range(5):
        for c in range(5):
            print(mat_numpy[l, c], end=" ")
        print()

    error = False
    for l in range(5):
        for c in range(5):
            if mat_numpy[l, c] != mat_om.value(l, c):
                print("matrices differ at:", l, c)
                error = True

    if error:
        print("conversion between OpenMEEG:Matrix <> numpy.ndarray is OK")

    # Testing going back from OpenMEEG to numpy
    mat_om2np = mat_om.array()
    np.testing.assert_array_equal(mat_numpy, mat_om2np)
