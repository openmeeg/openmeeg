import os.path as op

import numpy as np
import pytest

import openmeeg as om


def test_sensors_bad_matrix_layout_raises():
    """A non-Fortran-order Matrix argument must raise, not crash (gh-584)."""
    labels = ["toto", "tata"]
    positions = np.array([[0.0, 1.0, 2.0], [0.0, 1.0, 2.0]])  # C order, not F
    assert not positions.flags["F_CONTIGUOUS"]
    orientations = np.array([[-1.0, -1.0, -2.0], [-1.0, -1.0, -2.0]], order="F")
    weights = np.array([0.5, 0.5])
    radii = np.array([1.0, 1.0])
    with pytest.raises(Exception, match="Fortran order"):
        om.Sensors(labels, positions, orientations, weights, radii)


def test_sensors_non_contiguous_vector():
    """Read non-contiguous weights/radii correctly, not silently misread (gh-584)."""
    labels = ["toto", "tata"]
    positions = np.array([[0, 1, 2], [0, 1, 2]], order="F", dtype=float)
    orientations = np.array([[-1, -1, -2], [-1, -1, -2]], order="F", dtype=float)
    # Slice out of larger arrays so weights/radii are non-contiguous views.
    weights = np.array([[0.5, -1.0], [0.75, -1.0]])[:, 0]
    radii = np.array([[1.0, -1.0], [2.0, -1.0]])[:, 0]
    assert not weights.flags["C_CONTIGUOUS"]
    assert not radii.flags["C_CONTIGUOUS"]
    s = om.Sensors(labels, positions, orientations, weights, radii)
    np.testing.assert_array_equal(s.getWeights().array(), weights)
    np.testing.assert_array_equal(s.getRadii().array(), radii)


def test_sensors(data_path):
    # Test creating sensors from numpy arrays
    labels = ["toto", "tata"]
    positions = np.array([[0, 1, 2], [0, 1, 2]], order="F")
    orientations = np.array([[-1, -1, -2], [-1, -1, -2]], order="F")
    weights = np.array([0.5, 0.5])
    radii = np.array([1, 1])
    s1 = om.Sensors(labels, positions, orientations, weights, radii)
    print("s1 =", s1)
    s1.info()

    # Test creating sensors from only positions and orientations
    # XXX it's still broken for the moment
    # s1 = om.Sensors(positions, orientations)

    # Test creating Sensors from a geometry
    subject = "Head1"
    dirpath = op.join(data_path, subject)
    geom = om.read_geometry(
        op.join(dirpath, subject + ".geom"), op.join(dirpath, subject + ".cond")
    )
    om.Sensors(geom)

    # Test creating Sensors from file and a geometry
    om.Sensors(op.join(dirpath, subject + ".patches"), geom)

    # Test creating Sensors from an array and a geometry
    # XXX should be possible to pass positions as array directly
    om.Sensors(om.Matrix(positions), geom)
