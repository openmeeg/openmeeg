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

    # Test creating MEG sensors from only positions and orientations
    s2 = om.make_sensors(positions, orientations)
    assert s2.getNumberOfSensors() == 2

    # Test creating Sensors from a geometry
    subject = "Head1"
    dirpath = op.join(data_path, subject)
    geom = om.read_geometry(
        op.join(dirpath, subject + ".geom"), op.join(dirpath, subject + ".cond")
    )
    om.Sensors(geom)

    # Test creating Sensors from file and a geometry
    om.Sensors(op.join(dirpath, subject + ".patches"), geom)

    # Test creating EEG sensors from an array and a geometry (no om.Matrix wrap)
    om.make_sensors(positions, geometry=geom)


def test_make_sensors(data_path):
    """The make_sensors helper builds MEG/EEG sensors from plain arrays."""
    rng = np.random.RandomState(0)
    n = 5
    # C-ordered arrays (make_sensors converts to Fortran internally)
    positions = rng.randn(n, 3) * 0.02 + [0, 0, 0.12]
    orientations = positions / np.linalg.norm(positions, axis=1, keepdims=True)
    assert not positions.flags["F_CONTIGUOUS"]

    # MEG from positions + orientations, with defaults
    sensors = om.make_sensors(positions, orientations)
    assert sensors.getNumberOfSensors() == n
    np.testing.assert_allclose(sensors.getPositions().array(), positions)
    np.testing.assert_allclose(sensors.getWeights().array(), 1.0)
    np.testing.assert_allclose(sensors.getRadii().array(), 0.0)

    # Custom labels / weights / radii
    sensors = om.make_sensors(
        positions,
        orientations,
        labels=[f"MEG{i}" for i in range(n)],
        weights=np.full(n, 2.0),
        radii=np.full(n, 0.01),
    )
    np.testing.assert_allclose(sensors.getWeights().array(), 2.0)

    # EEG electrodes from positions + geometry
    subject = "Head1"
    dirpath = op.join(data_path, subject)
    geom = om.read_geometry(
        op.join(dirpath, subject + ".geom"), op.join(dirpath, subject + ".cond")
    )
    eeg = om.make_sensors(positions, geometry=geom)
    assert eeg.getNumberOfSensors() == n

    # Bad shapes raise cleanly (gh-584: previously segfaulted)
    with pytest.raises(ValueError, match="shape .n_sensors, 3."):
        om.make_sensors(np.zeros((n, 2)), orientations)
    with pytest.raises(ValueError, match="same shape as positions"):
        om.make_sensors(positions, np.zeros((n, 2)))
    with pytest.raises(ValueError, match=f"shape .{n},."):
        om.make_sensors(positions, orientations, weights=np.ones(n + 1))
    with pytest.raises(ValueError, match=f"length {n}"):
        om.make_sensors(positions, orientations, labels=["a", "b"])
    with pytest.raises(ValueError, match="require a geometry"):
        om.make_sensors(positions)
