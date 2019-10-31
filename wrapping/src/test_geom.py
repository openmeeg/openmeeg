#!/usr/bin/env python

import sys,os
import openmeeg as om
import numpy as np

Vertices  = np.array([[0.0,0.0,0.0],[1.0,0.0,0.0],
                     [1.0,1.0,0.0],[0.0,1.0,0.0]])
Triangles = np.array([[1,2,3],[2,3,0]])

mesh = om.Mesh(Vertices,Triangles)

g = om.Geometry( )

assert g.check(mesh) == True

g.import_meshes( [ mesh ] )

assert g.check(mesh) == False
