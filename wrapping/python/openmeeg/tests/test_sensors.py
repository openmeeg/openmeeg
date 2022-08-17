import numpy as np
import openmeeg as om


def test_sensors():
    labels = ["toto"]
    positions = np.array([[0, 1, 2], [0, 1, 2]], order="F")
    orientations = np.array([[-1, -1, -2], [-1, -1, -2]], order="F")
    weights = np.array([0.5, 0.5])
    radii = np.array([1, 1])
    s1 = om.Sensors(labels, positions, orientations, weights, radii)
    print("s1 =", s1)
    s1.info()
