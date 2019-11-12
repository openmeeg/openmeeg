#!/usr/bin/env python

import sys,os
import openmeeg as om
import numpy as np

data_path = path.dirname(path.abspath(__file__))
parser = OptionParser()
parser.add_option("-p","--path",dest="data_path",help="path to data folder",metavar="FILE",default=data_path)

options,args = parser.parse_args()
data_path    = options.data_path

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

BadVertices = np.array([0.0,0.0])
test_mesh("6",BadVertices,Triangles,False)

BadTriangles = np.array([1,2,3])
test_mesh("7",Vertices,BadTriangles,False)

Triangles = np.array([[1,2,3],[2,3,0]], dtype = np.uint64)
test_mesh("8",Vertices,Triangles,True)

Triangles = np.array([[1,2,3],[2,3,0]], dtype = np.int64)
test_mesh("9",Vertices,Triangles,True)

# test X -> should be OK
# TODO: Does not work if not jls....
data_file = os.path.join(data_path,"Head1","Head1.tri")
mesh_X = om.Mesh()
mesh_X.load(data_file)
mesh_X.update()

# test Y -> redo with np.array()
V_Y = mesh_X.vertices()
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
