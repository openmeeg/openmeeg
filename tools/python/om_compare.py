# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 2013

@author: - E. Olivi
"""

import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
from numpy import linalg as la

def rdmmag(X1,X2):
    X1 = vtk_to_numpy(X1)
    X2 = vtk_to_numpy(X2)
    rdm  = la.norm(X1/la.norm(X1)-X2/la.norm(X2))
    rmag = abs(1.-la.norm(X2)/la.norm(X1))
    print "  RDM = ",rdm, "\t  rMAG = ",rmag
    return rdm, rmag

def om_compare_vtp(f1,f2):
    """
    This function computes the rdm and (relative)mag errors of VTK::vtp files generated with OpenMEEG.
    Such a file defines a polydata, containing points and triangles of several
    meshes which are labelled through a vtkAbstractArray (Strings) associated to
    the cells (mesh names).
    Results of the forward problem (or a cortical mapping) can be seen thanks to
    arrays associated to points and cells (respectively potentials and normals
    currents).
    """
    maxrdm = 4
    maxmag = 2e500
    # Read the files
    reader1 = vtk.vtkXMLPolyDataReader()
    reader2 = vtk.vtkXMLPolyDataReader()
    reader1.SetFileName(f1); reader1.Update()
    poly1 = reader1.GetOutput()
    reader2.SetFileName(f2); reader2.Update()
    poly2 = reader2.GetOutput()
    # determine number of sources
    nb_sources = poly1.GetPointData().GetNumberOfArrays()-1
    # ensure it is the same geometry
    assert(poly1.GetNumberOfPoints() == poly2.GetNumberOfPoints())
    assert(poly1.GetNumberOfCells() == poly2.GetNumberOfCells())
    # Get the mesh names
    cell_labels = poly1.GetCellData().GetAbstractArray(0)
    assert(cell_labels.GetName() == 'Names')
    s = set(); nb_meshes = 0; cell_ids = list()
    for i in range(cell_labels.GetNumberOfValues()):
        s.add(cell_labels.GetValue(i))
        if len(s)>nb_meshes:
            # if a label is added, store the ID for connectivity filter
            cell_ids.append(i)
            nb_meshes += 1

    # TODO remove the min here for full comparisons
    for i_s in range(min(nb_sources,5)):
        print'\033[91m' + "Source " + str(i_s) + " :" + '\033[0m'
        V1 = poly1.GetPointData().GetArray('Potentials-'+str(i_s));
        P1 = poly1.GetCellData().GetArray('Currents-'+str(i_s));
        V2 = poly2.GetPointData().GetArray('Potentials-'+str(i_s));
        P2 = poly2.GetCellData().GetArray('Currents-'+str(i_s));
        print '\033[94m'+" all potentials :",
        rd0, mg0 = rdmmag(V1,V2)
        print " all currents   :",
        rd1, mg1 = rdmmag(P1,P2)
        maxrdm = min(maxrdm, rd0+rd1)
        maxmag = min(maxmag, mg0+mg1)
        print '\033[0m',
        for i_m in range(nb_meshes):
            print(" On Mesh " + cell_labels.GetValue(cell_ids[i_m]) + " :")
            conn1 = vtk.vtkPolyDataConnectivityFilter(); conn1.SetInput(poly1)
            conn1.SetExtractionModeToCellSeededRegions(); conn1.AddSeed(cell_ids[i_m]); conn1.Update()
            conn2 = vtk.vtkPolyDataConnectivityFilter(); conn2.SetInput(poly2)
            conn2.SetExtractionModeToCellSeededRegions(); conn2.AddSeed(cell_ids[i_m]); conn2.Update()
            V1 = conn1.GetOutput().GetPointData().GetArray('Potentials-'+str(i_s));
            P1 = conn1.GetOutput().GetCellData().GetArray('Currents-'+str(i_s));
            V2 = conn2.GetOutput().GetPointData().GetArray('Potentials-'+str(i_s));
            P2 = conn2.GetOutput().GetCellData().GetArray('Currents-'+str(i_s));
            print " \t potentials :",
            rdmmag(V1,V2)
            if (cell_labels.GetValue(cell_ids[i_m]) != "Scalp")&(cell_labels.GetValue(cell_ids[i_m]) != "scalp"):
                print " \t currents   :",
                rd,mg = rdmmag(P1,P2)
    return maxrdm, maxmag
