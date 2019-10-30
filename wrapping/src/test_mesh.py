#!/usr/bin/env python

import sys,os
import openmeeg as om
import numpy as np

def test_mesh(name,Vs,Ts,ExpectedResult):
    try:
        mesh = om.Mesh(Vs,Ts)
        mesh.update()
        mesh.info()
    except:
        if ExpectedResult:
            print("Test",name,"--> Failed")
            assert False
        else:
            print("Test",name,"--> Expected failure")
        return
    if not ExpectedResult:
        print("Test",name,"--> Unexpected success")
        assert False
    print("Test",name,"--> Expected success")

Vertices  = np.array([[0.0,0.0,0.0],[1.0,0.0,0.0],
                     [1.0,1.0,0.0],[0.0,1.0,0.0]])
Triangles = np.array([[1,2,3],[2,3,0]])

BadVertices = np.array([[0.0,0.0],    [1.0,0.0,0.0],
                        [1.0,1.0,0.0],[0.0,1.0,0.0]])
BadTriangles = np.array([[1,2  ],[0,1,2]])

test_mesh("1",Vertices,Triangles,True);
test_mesh("2",np.array([0.0,0.0,0.0]),Triangles,False)
test_mesh("3",Vertices,np.array([0,1,2]),False)
test_mesh("4",BadVertices,Triangles,False)
test_mesh("5",Vertices,BadTriangles,False)

# test 6 -> should be OK
# TODO: Does not work if not jls....
test_dir  = os.path.dirname(os.path.abspath(__file__))
data_file = os.path.join(test_dir,"..","..","..","data/Head1/Head1.tri")
mesh_6 = om.Mesh()
mesh_6.load(data_file)

# test 7 -> redo with np.array()
V6 = mesh_6.vertices()
#T6 = mesh_6.triangles()
#mesh_7 = om.Mesh(V6, T6)
#mesh_7.info()

# TODO
#
# mesh_6.nb_vertices()  == mesh_7.nb_vertices()
# mesh_6.nb_triangles() == mesh_7.nb_triangles()
# V7 = mesh_6.vertices()
# T7 = mesh_6.triangles()
# V6 == V7
# T6 == T7
