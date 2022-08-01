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

namespace OpenMEEG::GeometryIOs {

    /// \brief Handling of VTK\\vtp files.

#ifdef USE_VTK

    class OPENMEEG_EXPORT Vtp: public GeometryIO {

        typedef GeometryIO base;

    public:

        /// \brief load a VTK\\vtp file \param filename into a mesh. Optionally read some associated data in matrix \param data if \param READ_DATA is true.

        void load_meshes(Geometry& geometry) override {
            vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
            reader->SetFileName(fname.c_str());
            reader->Update();

            vtkMesh = reader->GetOutput();

            const unsigned npts  = vtkMesh->GetNumberOfPoints();
            const unsigned ntrgs = vtkMesh->GetNumberOfCells();

            //  Check that the mesh is a triangulation.

            for (unsigned i=0; i<ntrgs; ++i)
                if (vtkMesh->GetCellType(i)!=VTK_TRIANGLE) {
                    std::cerr << "This is not a triangulation" << std::endl;
                    exit(1);
                }
            
            //  Add vertices.

            int trash;
            vtkSmartPointer<vtkUnsignedIntArray> v_indices = vtkUnsignedIntArray::SafeDownCast(vtkMesh->GetPointData()->GetArray("Indices",trash));

            geometry.vertices().reserve(npts);

            IndexMap indmap;
            for (unsigned i=0; i<npts; ++i) {
                Vertex vertex(vtkMesh->GetPoint(i)[0],vtkMesh->GetPoint(i)[1],vtkMesh->GetPoint(i)[2]);
                if (trash!=-1) // index provided
                    vertex.index() = v_indices->GetValue(i);
                indmap.insert({ i, geometry.add_vertex(vertex) });
            }

            // Create meshes using the different meshes names associated with the cells.

            vtkSmartPointer<vtkStringArray> cell_id = vtkStringArray::SafeDownCast(vtkMesh->GetCellData()->GetAbstractArray("Names"));

            std::set<std::string> mesh_names;
            for (unsigned i=0; i<ntrgs; ++i) {
                const std::string& name = cell_id->GetValue(i);
                if (mesh_names.insert(name).second)
                    geometry.add_mesh(name);
            }

            // Insert the triangle and mesh vertices addresses into the right mesh

            int indices;
            vtkSmartPointer<vtkUnsignedIntArray> c_indices = vtkUnsignedIntArray::SafeDownCast(vtkMesh->GetCellData()->GetArray("Indices",indices));
            for (unsigned i=0; i<ntrgs; ++i) {
                const std::string& mesh_name = cell_id->GetValue(i);
                Mesh& mesh = geometry.mesh(mesh_name);
                const vtkSmartPointer<vtkIdList>& vl = vtkMesh->GetCell(i)->GetPointIds();
                const TriangleIndices t(vl->GetId(0),vl->GetId(1),vl->GetId(2));
                Triangle& triangle = mesh.add_triangle(t);
                if (indices!=-1) // index provided
                    triangle.index() = c_indices->GetValue(i);
            }

            for (auto& mesh : geometry.meshes()) {
                // Sets do not preserve the order, and we would like to preserve it so we push_back in the vector as soon as the element is unique.

                std::set<const Vertex*> mesh_vertices;
                mesh.vertices().clear();
                for (auto& triangle : mesh.triangles())
                    for (auto& vertex : triangle)
                        if (mesh_vertices.insert(vertex).second)
                            mesh.vertices().push_back(vertex);

                mesh.update(true);
            }
        }

        Matrix load_data() const override {

            const unsigned nl = vtkMesh->GetNumberOfPoints()+vtkMesh->GetNumberOfCells(); 
            const unsigned nc = vtkMesh->GetPointData()->GetNumberOfArrays()-1;

            if (nl==0 || nc==0)
                return Matrix();

            Matrix data(nl,nc);
            for (unsigned j=0; j<nc; ++j) {
                std::ostringstream sic;
                sic << j;
                int trash;
                const std::string& pot_array_name = "Potentials-"+sic.str();
                const unsigned potsz = vtkMesh->GetPointData()->GetArray(pot_array_name.c_str(),trash)->GetSize();
                for (unsigned i=0; i<potsz; ++i) {
                    const unsigned ind = vtkUnsignedIntArray::SafeDownCast(vtkMesh->GetPointData()->GetArray("Indices",trash))->GetValue(i);
                    data(ind,j) = vtkDoubleArray::SafeDownCast(vtkMesh->GetPointData()->GetArray(pot_array_name.c_str(),trash))->GetValue(i);
                }
                const std::string& curr_array_name = "Currents-"+sic.str();
                const unsigned currsz = vtkMesh->GetCellData()->GetArray(curr_array_name.c_str(),trash)->GetSize();
                for (unsigned i=0; i<currsz; ++i) {
                    const unsigned ind = vtkUnsignedIntArray::SafeDownCast(vtkMesh->GetCellData()->GetArray("Indices",trash))->GetValue(i);
                    data(ind,j) = vtkDoubleArray::SafeDownCast(vtkMesh->GetCellData()->GetArray(curr_array_name.c_str(),trash))->GetValue(i);
                }
            }
            return data;
        }

        /// \brief write a VTK\\vtp file.

        void save_geom(const Geometry& geometry) override {

            vtkMesh = vtkSmartPointer<vtkPolyData>::New();

            //  Vertices

            std::map<const Vertex*,unsigned> map;
            vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkUnsignedIntArray> point_indices = vtkSmartPointer<vtkUnsignedIntArray>::New();
            point_indices->SetName("Indices");

            unsigned i = 0;
            for (const auto& vertex : geometry.vertices()) {
                points->InsertNextPoint(vertex(0),vertex(1),vertex(2));
                point_indices->InsertNextValue(vertex.index());
                map[&vertex] = i++;
            }

            vtkMesh->SetPoints(points);
            vtkMesh->GetPointData()->AddArray(point_indices);

            //  Normals.

            vtkSmartPointer<vtkDoubleArray> normals = vtkSmartPointer<vtkDoubleArray>::New();
            normals->SetNumberOfComponents(3); // 3d normals (ie x, y, z)
            normals->SetName("Normals");

            // TODO: Not finished.

            //  Triangles.

            vtkSmartPointer<vtkCellArray> cells  = vtkSmartPointer<vtkCellArray>::New(); // the triangles

            vtkSmartPointer<vtkStringArray> cell_id = vtkSmartPointer<vtkStringArray>::New(); // ids/mesh names
            cell_id->SetName("Names");

            vtkSmartPointer<vtkUnsignedIntArray> cell_indices = vtkSmartPointer<vtkUnsignedIntArray>::New(); // indices
            cell_indices->SetName("Indices");

            for(const auto& mesh : geometry.meshes())
                for (const auto& triangle: mesh.triangles()) {
                    const vtkIdType vtktriangle[3] = { map[&(triangle.vertex(0))], map[&(triangle.vertex(1))], map[&(triangle.vertex(2))] };
                    cells->InsertNextCell(3,vtktriangle);
                    cell_id->InsertNextValue(mesh.name());
                    cell_indices->InsertNextValue(triangle.index());
                }

            vtkMesh->SetPolys(cells);
            vtkMesh->GetCellData()->AddArray(cell_id);
            vtkMesh->GetCellData()->AddArray(cell_indices);

        }

        void save_data(const Geometry& geometry,const Matrix& data) const override {
            // Check the data corresponds to the geometry

            const bool outer_interface_used = data.nlin()!=geometry.meshes().size();
            if (outer_interface_used )
                om_error(data.nlin()==geometry.meshes().size()-geometry.outermost_interface().nb_triangles());

            std::vector<vtkSmartPointer<vtkDoubleArray>> potentials(data.ncol()); // potential on vertices
            std::vector<vtkSmartPointer<vtkDoubleArray>> currents(data.ncol()); // current on triangles

            for (unsigned j = 0; j < data.ncol(); ++j) {
                std::stringstream sdip;
                sdip << j;

                potentials[j] = vtkSmartPointer<vtkDoubleArray>::New();
                potentials[j]->SetName(("Potentials-"+sdip.str()).c_str());

                currents[j]   = vtkSmartPointer<vtkDoubleArray>::New();
                currents[j]->SetName(("Currents-"+sdip.str()).c_str());

                for (const auto& vertex : geometry.vertices())
                    potentials[j]->InsertNextValue(data(vertex.index(),j));

                for(const auto& mesh : geometry.meshes())
                    for (const auto& triangle: mesh.triangles())
                        currents[j]->InsertNextValue(((mesh.outermost() && !outer_interface_used ) ? 0.0 : data(triangle.index(),j)));
            }

            for (unsigned j=0;j<data.ncol();++j) {
                vtkMesh->GetPointData()->AddArray(potentials[j]);
                vtkMesh->GetCellData()->AddArray(currents[j]);
            }
        }

        void write() const override {
            vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
            writer->SetFileName(fname.c_str());

            #if VTK_MAJOR_VERSION <= 5
            writer->SetInput(vtkMesh);
            #else
            writer->SetInputData(vtkMesh);
            #endif
            writer->Write();
        }

    private:

        const char* name() const override { return "vtp"; }

        GeometryIO* clone(const std::string& filename) const override { return new Vtp(filename); }

        Vtp(const std::string& filename=""): base(filename,"vtp") { }

        static const Vtp prototype;

        vtkSmartPointer<vtkPolyData> vtkMesh;
    };
#else

    class OPENMEEG_EXPORT Vtp: public GeometryIO {

        typedef GeometryIO base;

    public:

        const char* name() const override { return "vtp"; }

    protected:

        void load_meshes(Geometry&) override {
            std::cerr << "Error: please specify USE_VTK to cmake" << std::endl; // TODO in Legacy format ? // Exceptions
        }

        Matrix load_data() const override {
            std::cerr << "Error: please specify USE_VTK to cmake" << std::endl; // TODO in Legacy format ? // Exceptions
            throw "";
        }

    private:

        virtual void save_geom(const Geometry& geometry)                            override { }
        virtual void save_data(const Geometry& geometry,const Matrix& matrix) const override { }
        virtual void write()                                                  const override { }

        GeometryIO* clone(const std::string& filename) const override { return new Vtp(filename); }

        Vtp(const std::string& filename=""): base(filename,"vtp") { }

        static const Vtp prototype;
    };
#endif
}
