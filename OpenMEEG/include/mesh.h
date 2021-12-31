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

#include <vector>
#include <map>
#include <string>
#include <memory>

#include <om_common.h>
#include <triangle.h>
#include <om_utils.h>

#include <symmatrix.h>
#include <block_matrix.h>

namespace OpenMEEG {

    class Geometry;
    using maths::Range;
    using maths::Ranges;

    //  Mesh class
    //  \brief Mesh is a collection of triangles associated to a geometry containing the points
    //  on which triangles are based.

    class OPENMEEG_EXPORT Mesh {

        static Geometry* create_geometry(Geometry* geom);

    public:

        friend class Geometry;
        friend class MeshIO;

        typedef std::map<const Vertex*,TrianglesRefs> VertexTriangles;

        /// Default constructor
        /// or constructor using a provided geometry \param geometry

        Mesh(Geometry* geometry=nullptr): geom(create_geometry(geometry)) { }

        /// Constructor from scratch (vertices/triangles t be added)
        /// \param nv space to allocate for vertices
        /// \param nt space to allocate for triangles
        /// \param geometry the geometry to use
        // Do we need this ? TODO

        Mesh(const unsigned nv,const unsigned nt,Geometry* geometry=nullptr);

        #ifndef WIN32
        Mesh(const Mesh&) = delete;
        Mesh(Mesh&& m) = default;
        #endif

        /// Constructors
        /// \param filename mesh file name
        /// \param verbose verbose mode
        /// \param geometry geometry

        Mesh(const std::string& filename,const bool verbose,Geometry* geometry=nullptr): Mesh(geometry) {
            load(filename,verbose);
        }

        /// \param filename mesh file name
        /// \param geometry geometry

        Mesh(const std::string& filename,Geometry* geometry=nullptr): Mesh(filename,false,geometry) { }

        /// Destructor

        ~Mesh() { clear(); }

        std::string&       name()       { return mesh_name; } ///< \return the mesh name
        const std::string& name() const { return mesh_name; } ///< \return the mesh name

              VerticesRefs& vertices()       { return mesh_vertices; } ///< \return the vector of pointers to the mesh vertices
        const VerticesRefs& vertices() const { return mesh_vertices; } ///< \return the vector of pointers to the mesh vertices

              Geometry&    geometry() const { return *geom; }

              Triangles& triangles()       { return mesh_triangles; } ///< \return the triangles of the mesh
        const Triangles& triangles() const { return mesh_triangles; } ///< \return the triangles of the mesh

        TriangleIndices triangle(const Triangle& t) const;

        bool  current_barrier() const { return current_barrier_; }
        bool& current_barrier()       { return current_barrier_; }
        bool  isolated()        const { return isolated_;        }
        bool& isolated()              { return isolated_;        }

        /// \brief Add a triangle specified by its indices in the geometry.

        Triangle& add_triangle(const TriangleIndices inds);

        Triangle& add_triangle(const TriangleIndices inds,const IndexMap& indmap) {
            const TriangleIndices t = { indmap.at(inds[0]), indmap.at(inds[1]), indmap.at(inds[2])};
            return add_triangle(t);
        }

        void add(const std::vector<TriangleIndices>& trgs) {
            for (const auto& triangle : trgs)
                add_triangle(triangle);
        }

        void add(const std::vector<TriangleIndices>& trgs,const IndexMap& indmap) {
            for (const auto& triangle : trgs)
                add_triangle(triangle,indmap);
        }

        bool operator==(const Mesh& m) const { return triangles()==m.triangles(); }
        bool operator!=(const Mesh& m) const { return triangles()!=m.triangles(); }

        /// \brief Print info
        ///  Print to std::cout some info about the mesh.
        ///  \return void \sa

        void info(const bool verbose=false) const; ///< \brief Print mesh information.
        bool has_self_intersection() const;        ///< \brief Check whether the mesh self-intersects.
        bool intersection(const Mesh&) const;      ///< \brief Check whether the mesh intersects another mesh.
        bool has_correct_orientation() const;      ///< \brief Check local orientation of mesh triangles.
        void generate_indices();                   ///< \brief Generate indices (if allocate).
        void update(const bool topology_changed);  ///< \brief Recompute triangles normals, area, and vertex triangles.
        void merge(const Mesh&,const Mesh&);       ///< Merge two meshes.

        /// \brief Get the ranges of the specific mesh in the global matrix.
        /// \return vector of Range \sa

        Ranges vertices_ranges() const {
            std::vector<size_t> indices;
            for (const auto& vertex : vertices())
                indices.push_back(vertex->index());
            std::sort(indices.begin(),indices.end());
            Ranges result;
            for (auto it=indices.begin(); it!=indices.end();) {
                auto it1 = it;
                for (auto it2=it1+1; it2!=indices.end() && *it2==*it1+1; it1=it2++);
                result.push_back(Range(*it,*it1));
                it = it1+1;
            }
            return result;
        }

        //  Triangles always have a contiguous range as they are never shared between meshes.

        Range triangles_range() const { return Range(triangles().front().index(),triangles().back().index()); }
        
        /// \brief Get the triangles adjacent to vertex \param V .

        TrianglesRefs triangles(const Vertex& V) const { return vertex_triangles.at(&V); }

        /// \brief Get the triangles adjacent to \param triangle .

        TrianglesRefs adjacent_triangles(const Triangle& triangle) const {
            std::map<Triangle*,unsigned> mapt;
            TrianglesRefs result;
            for (auto& vertex : triangle)
                for (const auto& t2 : triangles(*vertex))
                    if (++mapt[t2]==2)
                        result.push_back(t2);
            return result;
        }

        /// Change mesh orientation.

        void change_orientation() {
            for (auto& triangle : triangles())
                triangle.change_orientation();
        }

        void correct_local_orientation(); ///< \brief Correct the local orientation of the mesh triangles.
        void correct_global_orientation(); ///< \brief Correct the global orientation (if there is one).
        double solid_angle(const Vect3& p) const; ///< Given a point p, computes the solid angle of the mesh seen from \param p .
        Normal normal(const Vertex& v) const; ///< \brief Get normal at vertex.`
        void laplacian(SymMatrix &A) const; ///< \brief Compute mesh laplacian.

        bool& outermost()       { return outermost_; } /// \brief Returns True if it is an outermost mesh.
        bool  outermost() const { return outermost_; }

        /// \brief Smooth Mesh
        /// \param smoothing_intensity
        /// \param niter
        /// \return void

        void smooth(const double& smoothing_intensity, const unsigned& niter);

        /// \brief Compute the square norm of the surfacic gradient

        void gradient_norm2(SymMatrix &A) const;

        /// Read mesh from file
        /// \param filename can be .vtk, .tri (ascii), .off, .bnd or .mesh.
        /// Be verbose if \param verbose is true. The difference between
        /// read and load is that read just reads the file and does not update
        /// the geometry. Read has to be used when multiple meshes are used in
        /// a geometry. load reads a mesh.

        void load(const std::string& filename,const bool verbose=true);

        /// Save mesh to file
        /// \param filename can be .vtk, .tri (ascii), .bnd, .off or .mesh

        void save(const std::string& filename) const ;

    #ifndef SWIGPYTHON
    private:
    #endif

        //  This private method must be accessible from swig.

        void reference_vertices(const IndexMap& indmap); ///< \brief Construct mesh vertices references,

    private:

        /// Map edges to an integer

        typedef std::map<std::pair<const Vertex*,const Vertex*>,int> EdgeMap;

        void clear();

        /// Add the mesh \param m to the current mesh. Assumes that geometry vertices are properly
        /// reserved (i.e. the vector is not resized while adding the mesh.

        void add_mesh(const Mesh& m);

        // regarding mesh orientation

        const EdgeMap compute_edge_map() const;
        bool  triangle_intersection(const Triangle&,const Triangle&) const;

        /// P1gradient: aux function to compute the surfacic gradient

        Vect3 P1gradient(const Vect3& p0,const Vect3& p1,const Vect3& p2) const { return crossprod(p1,p2)/det(p0,p1,p2); }

        /// P0gradient_norm2: aux function to compute the square norm of the surfacic gradient

        double P0gradient_norm2(const Triangle& t1,const Triangle& t2) const {
            return sqr(dotprod(t1.normal(),t2.normal()))/(t1.center()-t2.center()).norm2();
        }

        // Create the map that for each vertex gives the triangles containing it.

        void make_adjacencies() {
            vertex_triangles.clear();
            for (auto& triangle : triangles())
                for (const auto& vertex : triangle)
                    vertex_triangles[vertex].push_back(&triangle);
        }

        typedef std::shared_ptr<Geometry> Geom;

        std::string      mesh_name = "";     ///< Name of the mesh.
        VertexTriangles  vertex_triangles;   ///< links[&v] are the triangles that contain vertex v.
        Geometry*        geom;               ///< Pointer to the geometry containing the mesh.
        VerticesRefs     mesh_vertices;      ///< Vector of pointers to the mesh vertices.
        Triangles        mesh_triangles;     ///< Vector of triangles.
        bool             outermost_ = false; ///< Is it an outermost mesh ? (i.e does it touch the Air domain)

        /// Multiple 0 conductivity domains

        bool             current_barrier_ = false;
        bool             isolated_        = false;
    };

    typedef std::vector<Mesh> Meshes;
}
