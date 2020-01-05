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

#include <iostream>
#include <fstream>

#include <map>
#include <string>

#include <om_utils.h>
#include <MeshIO.h>

#ifdef USE_VTK
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDataReader.h>
#include <vtkCellArray.h>
#include <vtkCharArray.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#endif

namespace OpenMEEG::MeshIOs {

    /// \brief Mesh io for VTK file format.

    class OPENMEEG_EXPORT Vtk: public MeshIO {

        typedef MeshIO base;

    public:

    #ifdef USE_VTK

        ~Vtk() {
            if (buffer)
                delete[] buffer;
        }

        void load_points(Geometry& geom) override {
            //  TODO: Replace with newer filesystem calls.
            fs.seekg (0,ios::end);
            const unsigned length = fs.tellg();
            fs.seekg (0,ios::beg);

            // Read data as a block and pass it to vtk.

            buffer = new char[length];
            fs.read(buffer,length);

            vtkSmartPointer<vtkCharArray> buf = vtkSmartPointer<vtkCharArray>::New();
            buf->SetArray(buffer,length,1);

            //  Create the vtk reader.
            //  Make it read from the InputArray 'buf' instead of the default.

            vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
            reader->SetInputArray(buf);
            reader->SetReadFromInputString(1);
            if (!reader->IsFilePolyData()) {
                std::cerr << "Mesh \"" << fname << "\" is not a valid vtk poly data file" << std::endl;
                reader->Delete();
                exit(1);
            }

            reader->Update();
            vtkMesh = reader->GetOutput();

            const unsigned npts = vtkMesh->GetNumberOfPoints();
            Vertices vertices;
            for (unsigned i=0;i<npts;++i)
                vertices.push_back(Vertex(vtkMesh->GetPoint(i)));
            indmap = geom.add_vertices(vertices);
        }

        void load_triangles(OpenMEEG::Mesh& mesh) override {
            reference_vertices(mesh);

            const unsigned ntrgs = vtkMesh->GetNumberOfCells();
            mesh.triangles().reserve(ntrgs);

            for (unsigned i = 0; i < ntrgs; ++i) {
                if (vtkMesh->GetCellType(i)!=VTK_TRIANGLE) {
                    std::cerr << "Mesh \"" << fname << "\" is not a triangulation" << std::endl;
                    exit(1);
                }
                vtkIdList* l = vtkMesh->GetCell(i)->GetPointIds();
                mesh.add_triangle(TriangleIndices(l->GetId(0),l->GetId(1),l->GetId(2)),indmap);
            }
        }
    #else
        static unsigned vtk_error() {
            std::cerr << "OpenMEEG was not compiled with VTK support. Specify USE_VTK in cmake." << std::endl;
            exit(1);
            return 0;
        }

        void     load_points(Geometry&)          override { vtk_error();        }
        void     load_triangles(OpenMEEG::Mesh&) override { vtk_error();        }
        void     load(OpenMEEG::Mesh&)           override { vtk_error();        }
    #endif

        void save(const OpenMEEG::Mesh& mesh,std::ostream& os) const override {
            os << "# vtk DataFile Version 2.0" << std::endl;
            os << "Mesh file generated by OpenMEEG" << std::endl;
            os << "ASCII" << std::endl;
            os << "DATASET POLYDATA" << std::endl;
            os << "POINTS " << mesh.vertices().size() << " float" << std::endl;

            const VertexIndices& vertex_index(mesh);
            for (const auto& vertex : mesh.vertices())
                os << *vertex << std::endl;

            os << "POLYGONS " << mesh.triangles().size() << ' ' << mesh.triangles().size()*4 << std::endl;
            for (const auto& triangle : mesh.triangles())
                os << "3 " << vertex_index(triangle,0) << ' '
                           << vertex_index(triangle,1) << ' '
                           << vertex_index(triangle,2) << std::endl;

            os << "CELL_DATA " << mesh.triangles().size() << std::endl;
            os << "POINT_DATA " << mesh.vertices().size() << std::endl;
            os << "NORMALS normals float" << std::endl;
            for (const auto& vertex : mesh.vertices())
                os << mesh.normal(*vertex) << std::endl;
        }

        MeshIO* clone(const std::string& filename) const override { return new Vtk(filename); }

    private:

        Vtk(const std::string& filename=""): base(filename,"vtk") { }

        #ifdef USE_VTK
        char*                        buffer = nullptr;
        vtkSmartPointer<vtkPolyData> vtkMesh;
        #endif
        
        static const Vtk prototype;

        const char* name() const override { return "VTK"; }
    };
}
