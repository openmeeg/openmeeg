import os.path as op
import numpy as np
import openmeeg as om


def test_sensors(data_path):
    # Test creating sensors from numpy arrays
    labels = ["toto"]
    positions = np.array([[0, 1, 2], [0, 1, 2]], order="F")
    orientations = np.array([[-1, -1, -2], [-1, -1, -2]], order="F")
    weights = np.array([0.5, 0.5])
    radii = np.array([1, 1])
    s1 = om.Sensors(labels, positions, orientations, weights, radii)
    print("s1 =", s1)
    s1.info()

    # Test creating Sensors from a geometry
    subject = "Head1"
    dirpath = op.join(data_path, subject)
    geom = om.Geometry(
        op.join(dirpath, subject + ".geom"), op.join(dirpath, subject + ".cond")
    )
    om.Sensors(geom)

    # Test creating Sensors from file and a geometry
    om.Sensors(op.join(dirpath, subject + ".patches"), geom)

    # Test creating Sensors from an array and a geometry
    # XXX should be possible to pass positions as array directly
    om.Sensors(om.Matrix(positions), geom)
