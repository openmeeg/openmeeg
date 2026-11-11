// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <iterator>
#include <string>
#include <vector>
#include <set>

#include <om_common.h>
#include <vertex.h>
#include <triangle.h>
#include <interface.h>
#include <domain.h>
#include <matrix.h>

#include <OMExceptions.H>

namespace OpenMEEG {

    class Mesh;

    /// \brief Geometry contains the electrophysiological model
    /// Vertices, meshes and domains are stored in this geometry.

    class OPENMEEG_EXPORT Geometry {
    public:

        using MeshParts       = std::vector<MeshList>;
        using MeshDefinition  = std::pair<std::string,std::string>; // First and second strings are respectively the name and the path of the mesh.
        using MeshDefinitions = std::vector<MeshDefinition>;

        /// Constructors

        Geometry() { }

        Geometry(const unsigned n) {
            meshes().reserve(n);
        }

        Geometry(const std::string& geomFileName,const bool OLD_ORDERING=false) {
            load(geomFileName,OLD_ORDERING);
        }

        Geometry(const std::string& geomFileName,const std::string& condFileName,const bool OLD_ORDERING=false) {
            load(geomFileName,condFileName,OLD_ORDERING);
        }

        //  Absolutely necessary or wrong constructor is called because of conversion of char* to bool.

        Geometry(const char* geomFileName,const bool OLD_ORDERING=false): Geometry(std::string(geomFileName),OLD_ORDERING) { }
        Geometry(const char* geomFileName,const char* condFileName,const bool OLD_ORDERING=false):
            Geometry(std::string(geomFileName),std::string(condFileName),OLD_ORDERING) { }

        void info(const bool verbose=false) const; ///< \brief Print information on the geometry
        bool has_conductivities()           const {
            for (const auto& domain : domains())
                if (!domain.has_conductivity())
                    return false;
            return true;
        }

        bool selfCheck()                    const; ///< \brief the geometry meshes intersect each other
        bool check(const Mesh& m)           const; ///< \brief check if m intersect geometry meshes
        bool check_inner(const Matrix& m)   const; ///< \brief check if dipoles are outside of geometry meshes

        void check_geometry_is_nested();

        bool is_nested() const { return nested; }

        /// \brief Return the list of vertices involved in the geometry.

              Vertices& vertices()       { return geom_vertices; }
        const Vertices& vertices() const { return geom_vertices; }

        /// \brief Add a vertex \param V to the geometry and return the index of V in the vector of vertices.

        unsigned add_vertex(const Vertex& V) {
            // Insert the vertex in the set of vertices if it is not already in.

            const Vertices::iterator vit = std::find(vertices().begin(),vertices().end(),V);
            if (vit!=vertices().end())
                return vit-vertices().begin();

            vertices().push_back(V);
            return vertices().size()-1;
        }

        Mesh& add_mesh(const std::string& name="") {

            //  It is dangerous to store the returned mesh because the vector can be reallocated.
            //  Use mesh(name) after all meshes have been added....

            meshes().emplace_back(this);
            Mesh& mesh = meshes().back();
            mesh.name() = name;
            return mesh;
        }

        IndexMap add_vertices(const Vertices& vs) {
            IndexMap indmap;
            for (unsigned i=0; i<vs.size(); ++i)
                indmap.insert({ i, add_vertex(vs[i]) });
            return indmap;
        }

        /// \brief Return the list of meshes involved in the geometry.

              Meshes& meshes()       { return geom_meshes; }
        const Meshes& meshes() const { return geom_meshes; }

        const MeshList&  interacting_meshes() const { return non_isolated_meshes; }
        const MeshPairs& adjacent_meshes()    const { return meshpairs;           }

        /// \brief returns the Mesh called \param name .

        Mesh& mesh(const std::string& name);

        /// \brief  Return the list of domains.

              Domains& domains()       { return geom_domains; }
        const Domains& domains() const { return geom_domains; }

        /// \brief Get specific domains.

        const Domain& domain(const std::string& name) const; ///< \brief returns the Domain called \param name
        const Domain& domain(const Vect3& p)          const; ///< \brief returns the Domain containing the point p \param p a point

        /// \brief  Return the list of domains containing a mesh.

        DomainsRefs domains(const Mesh& m) const {
            DomainsRefs result;
            for (const auto& domain : domains())
                if (domain.contains(m))
                    result.push_back(&domain);
            return result;
        }

        bool is_external(const Domain& domain) const { return outer_domain==&domain; }

        const Interface& external_interface()  const; ///< \brief returns the external interface (only valid for nested geometries).
        const Interface& innermost_interface() const; ///< \brief returns the innermost interface (only valid for nested geometries).

        const Interface& interface(const std::string& name) const; ///< \brief returns the Interface called \param name

        /// \brief Return the conductivity jump across a mesh (i.e. between the 2 domains it separates).

        double conductivity_jump(const Mesh& m) const {
            const DomainsRefs& doms = domains(m);
            double res = 0.0;
            for (const auto& domainptr : doms)
                res += domainptr->conductivity()*domainptr->mesh_orientation(m);
            return res;
        }

        //  Calling this method read induces failures due do wrong conversions when read is passed with one or two arguments...

        void load(const std::string& filename,const bool OLD_ORDERING=false) {
            clear();
            read_geometry_file(filename);
            finalize(OLD_ORDERING);
        }

        void load(const std::string& geomFileName,const std::string& condFileName,const bool OLD_ORDERING=false) {
            clear();
            read_geometry_file(geomFileName);
            read_conductivity_file(condFileName);
            finalize(OLD_ORDERING);
        }

        void set_conductivity(const std::string& domain_name,const double cond_value) {
            for (auto& domain : domains())
                if (domain.name()==domain_name) {
                    domain.set_conductivity(cond_value);
                    return;
                }
            throw OpenMEEG::BadDomain(domain_name);
        }

        void import(const MeshDefinitions& meshes);

        void save(const std::string& filename) const;

        void finalize(const bool OLD_ORDERING=false) {
            // TODO: We should check the correct decomposition of the geometry into domains here.
            // In a correct decomposition, each interface is used exactly twice with two opposite
            // Search for the external interface/meshes and set boolean external on each such mesh.
            // An external interface is defined as delimiting a maximal outside domain with null
            // conductivity. It is supposed to be unique. Typically it is the most outward mesh
            // delimiting a domain with a non-zero conductivity

            if (has_conductivities())
                set_mesh_properties(); // Mark meshes that touch the domains of null conductivity.

            if (domains().size()!=0) {
                set_external_domain(external_domain());
                check_geometry_is_nested();
            }

            generate_indices(OLD_ORDERING);
            make_mesh_pairs();
            #ifdef DEBUG
            for (const auto& mesh : meshes())
                mesh.check_consistency("geometry finalize step");
            #endif
        }

        /// Handle multiple isolated domains

        size_t  num_current_barrier_triangles() const { return nb_current_barrier_triangles; }
        size_t& num_current_barrier_triangles()       { return nb_current_barrier_triangles; }

        size_t num_parameters() const { return nb_params;                                        } ///< \brief the total number of vertices + triangles
        size_t num_unknowns()   const { return num_parameters()-num_current_barrier_triangles(); } ///< \brief the total number of unknowns

        const MeshParts& isolated_parts() const { return independant_parts; }

        void set_mesh_properties();

        const Mesh& mesh(const std::string& id) const; //  Is this useful ?? TODO.

    private:

        void clear() {
            geom_vertices.clear();
            geom_meshes.clear();
            geom_domains.clear();
            nested = false;
            outer_domain = 0;
            nb_params = 0;
        }

        void read_geometry_file(const std::string& filename);
        void read_conductivity_file(const std::string& filename);

        /// Returns the external domain.

        Domain& external_domain();

        void set_external_domain(Domain& domain) {
            outer_domain = &domain;
            for (auto& boundary : domain.boundaries())
                boundary.interface().set_to_external();
        }

        void make_mesh_pairs();

        /// Members

        Vertices      geom_vertices;
        Meshes        geom_meshes;
        Domains       geom_domains;

        const Domain* outer_domain = 0;
        bool          nested       = false;
        size_t        nb_params    = 0;   // total number = nb of vertices + nb of triangles

        void generate_indices(const bool);

        // Two meshes can have 0, 1 or 2 domains in common.
        // - 0: they are not interacting.
        // - 1: they are on two different boundaries of the domain.
        // - 2: they are on the same boundary of the domain.

        DomainsRefs common_domains(const Mesh& m1,const Mesh& m2) const {
            const DomainsRefs& doms1 = domains(m1);
            const DomainsRefs& doms2 = domains(m2);
            DomainsRefs doms;
            std::set_intersection(doms1.begin(),doms1.end(),doms2.begin(),doms2.end(),std::back_inserter(doms));
            return doms;
        }

        /// Handle multiple isolated domains.

        std::set<Vertex> invalid_vertices;                 ///< \brief  does not equal to the vertices of invalid meshes because there are shared vertices
        size_t           nb_current_barrier_triangles = 0; ///< \brief number of triangles with 0 normal current. Including triangles of invalid meshes.

        MeshParts independant_parts;   ///< \brief Meshes that belong to different isolated groups.
        MeshList  non_isolated_meshes; ///< \brief List of non isolated sub-meshes of the head.
        MeshPairs meshpairs;           ///< \brief Pairs of different and non isolated meshes.
    };
}
