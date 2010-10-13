# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 15:43:44 2010

@author: - B. Burle & A. Gramfort
"""

import numpy as np

class Mesh():
    '''
    A class allowing to open .tri files, store the corresponding mesh
    data (vertices and triangles), and plot their respective surfaces

    To read a .tri file, call the function 'self.read_tri(filename)'

    Example:
    --------
    head = Mesh()
    head.read_tri("MyHeadFile.tri")

    To plot the surface of the corresponding object, call the function
    'self.plot(**kwarg)'. The kwarg options are the one from
    enthought.mayavi.mlab.triangular_mesh

    '''
    def __init__(self, fname=None):
        self.points = []
        self.faces = []
        self.normals = []
        if fname is not None:
            self.read_tri(fname)

    def read_tri(self, fname):
        assert(fname.endswith('.tri'))
        fid = file(fname,"r")
        # read the number of vertices
        npoints = int(fid.readline().split()[1])
        # fills the vertices arrays
        for _ in xrange(npoints):
            vals = map(float,fid.readline().split())
            self.points.append(vals[:3])
            self.normals.append(vals[3:])

        # Read the number of triangles
        n_faces = int(fid.readline().split()[1])
        # create the list of triangles
        for _ in xrange(n_faces):
            vals = map(int,fid.readline().split())
            self.faces.append(vals[:3])

        # Convert to numpy arrays
        self.points = np.asarray(self.points)
        self.normals = np.asarray(self.normals)
        self.faces = np.asarray(self.faces)

    def plot(self,**kwarg):
        '''Plot mesh with Mayavi
        '''
        from enthought.mayavi import mlab
        f = mlab.triangular_mesh(self.points[:,0],self.points[:,1],self.points[:,2],self.faces, **kwarg)
        return f

