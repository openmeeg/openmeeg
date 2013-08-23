/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

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

#ifndef OPENMEEG_GEOMETRY_IO_H
#define OPENMEEG_GEOMETRY_IO_H

#include <map>
#ifdef USE_VTK
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
// #include <vtkPolyDataWriter.h>
#include <vtkDataReader.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkProperty.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkAbstractArray.h>
#include <vtkStringArray.h>
#include <vtkSmartPointer.h>
#include <vtkCleanPolyData.h>
#endif

namespace OpenMEEG {

    /// \brief load a VTK\vtp file \param filename
    void Geometry::load_vtp(const std::string& filename)
    {
    #ifdef USE_VTK
        vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(filename.c_str()); // Specify file name of vtp data file to read
        reader->Update();

        vtkSmartPointer<vtkPolyData> vtkMesh = reader->GetOutput();

        unsigned npts, ntrgs;
        npts = vtkMesh->GetNumberOfPoints();
        vertices_.reserve(npts); // alocate memory for the vertices

        if ( vtkMesh->GetPointData()->GetNormals()->GetNumberOfTuples() != vtkMesh->GetNumberOfPoints() ) {
            std::cout << "Recomputes the vertices normals..." << std::endl;
            vtkSmartPointer<vtkPolyDataNormals> newNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
            newNormals->SetInput(vtkMesh);
            newNormals->Update();
            vtkMesh = newNormals->GetOutput();
        }

        vtkSmartPointer<vtkDataArray> normals = vtkMesh->GetPointData()->GetNormals();

        if ( npts != normals->GetNumberOfTuples() ) {
            std::cerr << "Error: number of vertices is not equal to number of normals in vtp file, correct or remove the normals." << std::endl;
            exit(1);
        }

        if ( normals->GetNumberOfComponents() != 3 ) {
            std::cerr << "Error: wrong number of components of normals in vtp file, correct or remove the normals." << std::endl;
            exit(1);
        }

        for ( unsigned i = 0; i < npts; ++i) {
            vertices_.push_back(Vertex(vtkMesh->GetPoint(i)[0], vtkMesh->GetPoint(i)[1], vtkMesh->GetPoint(i)[2],
                                       normals->GetTuple(i)[0], normals->GetTuple(i)[1], normals->GetTuple(i)[2]));
        }

        // find the number of different meshes reading the string array associated with the cells
        std::set<std::string> meshes_name;

        // ids/mesh name
        vtkSmartPointer<vtkStringArray> cell_id = vtkStringArray::SafeDownCast(vtkMesh->GetCellData()->GetAbstractArray("Names"));

        ntrgs = vtkMesh->GetNumberOfCells();
        for ( unsigned i = 0; i < ntrgs; ++i) {
            meshes_name.insert(cell_id->GetValue(i));
        }

        // create the meshes
        for ( std::set<std::string>::const_iterator sit = meshes_name.begin(); sit != meshes_name.end(); ++sit) {
            meshes_.push_back(Mesh(vertices_, *sit));
        }

        // insert the triangle and mesh vertices address into the right mesh
        vtkSmartPointer<vtkIdList>   l;

        for ( iterator mit = begin(); mit != end(); ++mit) {
            for ( unsigned i = 0; i < ntrgs; ++i) {
                // get the mesh which has this name
                if ( cell_id->GetValue(i) == mit->name() ) {
                    if ( vtkMesh->GetCellType(i) == VTK_TRIANGLE ) {
                        l = vtkMesh->GetCell(i)->GetPointIds();
                        mit->push_back(Triangle(vertices_[l->GetId(0)],
                                                vertices_[l->GetId(1)],
                                                vertices_[l->GetId(2)])); 
                    } else {
                        std::cerr << "This is not a triangulation" << std::endl;
                        exit(1);
                    }
                }
            }
            mit->build_mesh_vertices();
            mit->update();
        }
    #else
        std::cerr << "Error: please specify USE_VTK to cmake" << std::endl; // TODO in Legacy format ? // Exceptions
    #endif
    }
    
    /// \brief write a VTK\vtp file \param filename
    void Geometry::write_vtp(const std::string& filename) const
    {

    #ifdef USE_VTK
        vtkSmartPointer<vtkPolyData>    polydata = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints>      points   = vtkSmartPointer<vtkPoints>::New();      // vertices
        vtkSmartPointer<vtkDoubleArray> normals  = vtkSmartPointer<vtkDoubleArray>::New(); // normals
        vtkSmartPointer<vtkStringArray> cell_id  = vtkSmartPointer<vtkStringArray>::New(); // ids/mesh name
        
        normals->SetNumberOfComponents(3); // 3d normals (ie x,y,z)
        normals->SetName("Normals");
        cell_id->SetName("Names");

        std::map<const Vertex *, unsigned> map;

        unsigned i = 0;
        for ( Vertices::const_iterator vit = vertex_begin(); vit != vertex_end(); ++vit, ++i) {
            points->InsertNextPoint((*vit)(0), (*vit)(1), (*vit)(2));
            const double n[3] = { vit->normal()(0), vit->normal()(1), vit->normal()(2)};
            normals->InsertNextTupleValue(n);
            map[&*vit] = i;
        }
        // Add the vertices to polydata
        polydata->SetPoints(points);
        // Add the normals to the points in the polydata
        polydata->GetPointData()->SetNormals(normals);
        
        vtkSmartPointer<vtkCellArray> polys  = vtkSmartPointer<vtkCellArray>::New(); // the triangles
        
        for ( Meshes::const_iterator mit = meshes_.begin(); mit != meshes_.end(); ++mit) {
            for ( Mesh::const_iterator tit = mit->begin(); tit != mit->end(); ++tit) {
                vtkIdType triangle[3] = { map[&(tit->s1())], map[&(tit->s2())], map[&(tit->s3())] };
                polys->InsertNextCell(3, triangle);
                cell_id->InsertNextValue(mit->name());
            }
        }

        // Add the triangles to polydata
        polydata->SetPolys(polys);
        // Add the array of Names to polydata
        polydata->GetCellData()->AddArray(cell_id);
        
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(filename.c_str());

        // writer->SetCompressorTypeToZLib();
        // writer->SetCompressorTypeToNone();

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

#endif  //! OPENMEEG_GEOMETRY_IO_H
