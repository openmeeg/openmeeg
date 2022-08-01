import os
from os import path as op
import pytest
import numpy as np

import openmeeg as om


@pytest.mark.skipif(os.getenv('OPENMEEG_BAD_MKL') == os.getenv('OPENMEEG_BAD_MSVC') == '1',
                    reason='bad windows mkl')
def test_python2(data_path):
    subject = "Head1"
    file_name_skeleton = op.join(data_path, subject, subject)

    # dipole
    dipoles = om.Matrix()
    dipoles.load(file_name_skeleton + ".dip")
    D = dipoles.array()
    print("D is a", D.__class__)
    print(D)
    # Examples of basic linear algebra
    print("Determinant of D is equal to: ", np.linalg.det(D))

    # TODO: sensors.read() using numpy arrays
    # TODO: sensors == [ double ]
    sensors = om.Sensors()
    sensors.load(file_name_skeleton + ".squids")
    # TODO: D = asarray(sensors).copy...
    # TODO: sensors_1 = om.Matrix(D)

    # TODO: mesh.read() using numpy arrays
    # TODO: mesh == [ double ] , [ int ]
    # mesh = om.Mesh()
    # mesh.load( file_name_skeleton + ".tri")
    # TODO: V = [...]
    # TODO: I = [...]
    # TODO: mesh_1 = om.Mesh(V, I)

    # TODO: geom.read() using raw python
    # TODO: conductivity == { string => double}
    # TODO: geometry     == { [ mesh ] , { string => [ int ] } }
    # geom = om.Geometry()
    # geom.read( file_name_skeleton + ".geom" , file_name_skeleton + ".cond" )
