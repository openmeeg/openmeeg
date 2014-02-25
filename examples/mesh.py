# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 15:43:44 2010

@author: - B. Burle & A. Gramfort
"""

import numpy as np
import openmeeg as om

class Mesh():
    '''
    A class allowing to open meshes
    data (vertices and triangles), and plot their respective surfaces

    Example:
    --------
    head = Mesh("MyHeadFile.tri")

    To plot the surface of the corresponding object, call the function
    'self.plot(**kwarg)'. The kwarg options are the one from
    mayavi.mlab.triangular_mesh

    '''
    def __init__(self, fname=None):
        self.points = []
        self.faces = []
        if fname is not None:
            self.load(fname)

    def load(self, fname):
        m = om.Mesh(fname)
        self.points = np.zeros((m.nb_vertices(),3),dtype=float)
        self.faces  = np.zeros((m.nb_triangles(),3),dtype=int)
        vit = m.vertex_begin(); iit = 0
        while vit != m.vertex_end():
            self.points[iit,:] = [vit.value()(0), vit.value()(1), vit.value()(2)]
            iit += 1
            vit.incr()
        tit = m.begin(); iit = 0
        while tit != m.end():
            self.faces[iit,:] = [tit.value().s1().getindex(), tit.value().s2().getindex(), tit.value().s3().getindex()]
            iit += 1
            tit.incr()

    def plot(self,**kwarg):
        '''Plot mesh with Mayavi
        '''
        from mayavi import mlab
        f = mlab.triangular_mesh(self.points[:,0],self.points[:,1],self.points[:,2],self.faces, **kwarg)
        return f

