/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#pragma once

#include <map>
#include <sstream>

#ifdef USE_VTK
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkAbstractArray.h>
#include <vtkStringArray.h>
#include <vtkSmartPointer.h>
#endif

namespace OpenMEEG {

    /// \brief load a VTK\\vtp file \param filename into a mesh. Optionally read some associated data in matrix \param data if \param READ_DATA is true.

    void Geometry::load_vtp(const std::string& filename, Matrix& data, const bool READ_DATA) {
        #ifdef USE_VTK
        vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(filename.c_str()); // Specify file name of vtp data file to read
        reader->Update();

        vtkSmartPointer<vtkPolyData> vtkMesh = reader->GetOutput();

        int trash;
        vtkSmartPointer<vtkUnsignedIntArray> v_indices = vtkUnsignedIntArray::SafeDownCast(vtkMesh->GetPointData()->GetArray("Indices", trash));
        vtkSmartPointer<vtkUnsignedIntArray> c_indices = vtkUnsignedIntArray::SafeDownCast(vtkMesh->GetCellData()->GetArray("Indices", trash));

        const unsigned npts = vtkMesh->GetNumberOfPoints();
        vertices_.reserve(npts); // alocate memory for the vertices

        for (unsigned i = 0; i < npts; ++i) {
            if (trash == -1) { // no indices provided
                vertices_.push_back(Vertex(vtkMesh->GetPoint(i)[0], vtkMesh->GetPoint(i)[1], vtkMesh->GetPoint(i)[2]));
            } else {
                vertices_.push_back(Vertex(vtkMesh->GetPoint(i)[0], vtkMesh->GetPoint(i)[1], vtkMesh->GetPoint(i)[2], v_indices->GetValue(i)));
            }
        }

        // find the number of different meshes reading the string array associated with the cells
        std::vector<std::string> meshes_name;

        // ids/mesh name
        vtkSmartPointer<vtkStringArray> cell_id = vtkStringArray::SafeDownCast(vtkMesh->GetCellData()->GetAbstractArray("Names"));

        const unsigned ntrgs = vtkMesh->GetNumberOfCells();
        for (unsigned i = 0; i < ntrgs; ++i) {
            if (std::find(meshes_name.begin(), meshes_name.end(), cell_id->GetValue(i)) == meshes_name.end()) {
                meshes_name.push_back(cell_id->GetValue(i));
            }
        }

        // create the meshes
        for (std::vector<std::string>::const_iterator vit = meshes_name.begin(); vit != meshes_name.end(); ++vit) {
            meshes_.push_back(Mesh(vertices_, *vit));
        }

        // insert the triangle and mesh vertices address into the right mesh
        vtkSmartPointer<vtkIdList> l;

        for (iterator mit = begin(); mit != end(); ++mit) {
            for (unsigned i = 0; i < ntrgs; ++i) {
                // get the mesh which has this name
                if (cell_id->GetValue(i) == mit->name()) {
                    if (vtkMesh->GetCellType(i) == VTK_TRIANGLE) {
                        l = vtkMesh->GetCell(i)->GetPointIds();
                        if (trash != -1) { // no indices provided
                            mit->push_back(Triangle(vertices_[l->GetId(0)],
                                                    vertices_[l->GetId(1)],
                                                    vertices_[l->GetId(2)],
                                                    c_indices->GetValue(i))); 
                        } else {
                            mit->push_back(Triangle(vertices_[l->GetId(0)],
                                                    vertices_[l->GetId(1)],
                                                    vertices_[l->GetId(2)]));
                        }
                    } else {
                        std::cerr << "This is not a triangulation" << std::endl;
                        exit(1);
                    }
                }
            }

            mit->build_mesh_vertices();
            mit->update();
            if (READ_DATA) {
                const unsigned nc = vtkMesh->GetPointData()->GetNumberOfArrays()-1;
                const unsigned nl = vtkMesh->GetNumberOfPoints()+vtkMesh->GetNumberOfCells(); 
                if (nc != 0) {
                    data = Matrix(nl, nc);
                    for (unsigned ic = 0; ic < nc; ++ic) {
                        std::ostringstream sic;
                        sic << ic;
                        for (unsigned iil = 0; iil < vtkMesh->GetPointData()->GetArray(("Potentials-"+sic.str()).c_str(), trash)->GetSize(); ++iil) {
                            data(vtkUnsignedIntArray::SafeDownCast(vtkMesh->GetPointData()->GetArray("Indices", trash))->GetValue(iil), ic) =
                                vtkDoubleArray::SafeDownCast(vtkMesh->GetPointData()->GetArray(("Potentials-"+sic.str()).c_str(), trash))->GetValue(iil);
                        }
                        for (unsigned iil = 0; iil < vtkMesh->GetCellData()->GetArray(("Currents-"+sic.str()).c_str(), trash)->GetSize(); ++iil) {
                            data(vtkUnsignedIntArray::SafeDownCast(vtkMesh->GetCellData()->GetArray("Indices", trash))->GetValue(iil), ic) =
                                vtkDoubleArray::SafeDownCast(vtkMesh->GetCellData()->GetArray(("Currents-"+sic.str()).c_str(), trash))->GetValue(iil);
                        }
                    }
                }
            }
        }
        #else
        std::cerr << "Error: please specify USE_VTK to cmake" << std::endl; // TODO in Legacy format ? // Exceptions
        #endif
    }

    /// \brief write a VTK\\vtp file \param filename with associated data in \param data .

    void Geometry::write_vtp(const std::string& filename, const Matrix& data) const {

        #ifdef USE_VTK
        vtkSmartPointer<vtkPolyData>         polydata = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints>           points   = vtkSmartPointer<vtkPoints>::New();      // vertices
        vtkSmartPointer<vtkDoubleArray>      normals  = vtkSmartPointer<vtkDoubleArray>::New(); // normals
        vtkSmartPointer<vtkStringArray>      cell_id  = vtkSmartPointer<vtkStringArray>::New(); // ids/mesh name
        vtkSmartPointer<vtkUnsignedIntArray> cell_indices = vtkSmartPointer<vtkUnsignedIntArray>::New(); // indices
        vtkSmartPointer<vtkUnsignedIntArray> point_indices = vtkSmartPointer<vtkUnsignedIntArray>::New(); // indices
        std::vector<vtkSmartPointer<vtkDoubleArray> > potentials(data.ncol()); // potential on vertices
        std::vector<vtkSmartPointer<vtkDoubleArray> > currents(data.ncol()); // current on triangles

        normals->SetNumberOfComponents(3); // 3d normals (ie x, y, z)
        normals->SetName("Normals");
        cell_id->SetName("Names");
        cell_indices->SetName("Indices");
        point_indices->SetName("Indices");

        if (data.nlin() != 0) {
            // Check the data corresponds to the geometry
            bool HAS_OUTERMOST = false; // data has last p values ?
            if (data.nlin() == size()) {
                HAS_OUTERMOST = true;
            } else {
                om_error(data.nlin() == size()-outermost_interface().nb_triangles());
            }
            for (unsigned j = 0; j < data.ncol(); ++j) {
                std::stringstream sdip;
                sdip << j;
                potentials[j] = vtkSmartPointer<vtkDoubleArray>::New();
                currents[j]   = vtkSmartPointer<vtkDoubleArray>::New();
                potentials[j]->SetName(("Potentials-"+sdip.str()).c_str());
                currents[j]->SetName(("Currents-"+sdip.str()).c_str());

                for (Vertices::const_iterator vit = vertex_begin(); vit != vertex_end(); ++vit) {
                    potentials[j]->InsertNextValue(data(vit->index(), j));
                }

                for (Meshes::const_iterator mit = meshes_.begin(); mit != meshes_.end(); ++mit) {
                    for (Mesh::const_iterator tit = mit->begin(); tit != mit->end(); ++tit) {
                        currents[j]->InsertNextValue(((mit->outermost() && !HAS_OUTERMOST) ? 0.0 : data(tit->index(), j)));
                    }
                }
            }
        }

        std::map<const Vertex*, unsigned> map;

        unsigned i = 0;
        for (Vertices::const_iterator vit = vertex_begin(); vit != vertex_end(); ++vit, ++i) {
            points->InsertNextPoint((*vit)(0), (*vit)(1), (*vit)(2));
            point_indices->InsertNextValue(vit->index());
            map[&*vit] = i;
        }
        // Add the vertices to polydata
        polydata->SetPoints(points);

        vtkSmartPointer<vtkCellArray> polys  = vtkSmartPointer<vtkCellArray>::New(); // the triangles

        for (Meshes::const_iterator mit = meshes_.begin(); mit != meshes_.end(); ++mit) {
            for ( Mesh::const_iterator tit = mit->begin(); tit != mit->end(); ++tit) {
                vtkIdType triangle[3] = { map[&(tit->s1())], map[&(tit->s2())], map[&(tit->s3())] };
                polys->InsertNextCell(3, triangle);
                cell_id->InsertNextValue(mit->name());
                cell_indices->InsertNextValue(tit->index());
            }
        }

        // Add the triangles to polydata
        polydata->SetPolys(polys);
        // Add the array of Names to polydata
        polydata->GetCellData()->AddArray(cell_id);
        // Add the array of Indices to polydata cells
        polydata->GetCellData()->AddArray(cell_indices);
        // Add the array of Indices to polydata points
        polydata->GetPointData()->AddArray(point_indices);
        // Add optional potentials and currents
        if (data.nlin() != 0) {
            for (unsigned j = 0; j < data.ncol(); ++j) {
                polydata->GetPointData()->AddArray(potentials[j]);
                polydata->GetCellData()->AddArray(currents[j]);
            }
        }

        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(filename.c_str());

        #if VTK_MAJOR_VERSION <= 5
        writer->SetInput(polydata);
        #else
        writer->SetInputData(polydata);
        #endif
        writer->Write();

        #else
        std::cerr << "Error: please specify USE_VTK to cmake" << std::endl; // TODO write in Legacy format ? // Exceptions
        #endif
    }
}
