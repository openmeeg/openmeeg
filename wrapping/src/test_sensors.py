#!/usr/bin/env python
import sys
import os
import numpy as np

import openmeeg as om


# nes Sensors construtor
Labels = ["toto"]
Positions = np.array([[0, 1, 2]])
Orientations = np.array([[-1, -1, -2]])
Weights = np.array([0.5])
Radii = np.array([1])

s1 = om.Sensors(Labels, Positions, Orientations, Weights, Radii)
print("s1 =", s1)
s1.info()
