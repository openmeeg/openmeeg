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

#include <geometry.h>
#include <MeshIO.h>
#include <GeometryIO.h>
#include <PropertiesSpecialized.h>

namespace OpenMEEG {

    /// Search for the outermost domain.

    Domain&
    Geometry::outermost_domain() {
        // The outermost domain is defined as the domain which has no inside.
        // This domain should be unique. Otherwise the decomposition of the geometry
        // into domains is wrong. We assume here that it is correct.

        for (auto& domain : domains()) {
            bool outer = true;
            for (auto& boundary : domain.boundaries())
                if (boundary.inside()) {
                    outer = false;
                    break;
                }
            if (outer)
                return domain;
        }

        warning("Geometry::outermost_domain: Error outermost domain is not defined.");
        throw OpenMEEG::BadDomain("outermost");
    }

    const Interface&
    Geometry::innermost_interface() const {

        if (!is_nested()) {
            warning("Geometry::innermost_interface: Error innermost interface is only defined for nested geometries.");
            throw OpenMEEG::BadInterface("innermost");
        }

        for (const auto& domain : domains()) {
            bool inner = true;
            for (const auto& boundary : domain.boundaries())
                if (!boundary.inside()) {
                    inner = false;
                    break;
                }
            if (inner)
                return domain.boundaries().front().interface();
        }

        // Should never append as this function should only be called for nested geometries.

        warning("Geometry::innerermost_interface: Error innermost interface is not defined.");
        throw OpenMEEG::BadInterface("innermost");
    }

    const Interface&
    Geometry::outermost_interface() const {
        for (const auto& domain : domains())
            if (is_outermost(domain))
                return domain.boundaries().front().interface();

        // Should never append

        warning("Geometry::outermost_interface: Error outermost interface were not set.");
        throw OpenMEEG::BadInterface("outermost");
    }

    Mesh& Geometry::mesh(const std::string& id) {
        for (auto& mesh: meshes())
            if (mesh.name()==id)
                return mesh;

        // Should never happen

        warning(std::string("Geometry::mesh: Error mesh id/name not found: ") + id);
        throw OpenMEEG::BadInterface(id);
    }

    const Mesh& Geometry::mesh(const std::string& id) const {
        for (const auto& mesh: meshes())
            if (mesh.name()==id)
                return mesh;

        // Should never happen

        warning(std::string("Geometry::mesh: Error mesh id/name not found: ") + id);
        throw OpenMEEG::BadInterface(id);
    }

    void Geometry::info(const bool verbose) const {
        if (is_nested()) {
            std::cout << "This geometry is a NESTED geometry." << std::endl;
        } else {
            int shared = -vertices().size();
            for (const auto& mesh : meshes())
                shared += mesh.vertices().size();

            // those are not the number of shared vertices but the number of demands for adding the same vertex...
            std::cout << "This geometry is a NON NESTED geometry. (There was " << shared << " demands for adding same vertices)." << std::endl;
        }

        for (const auto& mesh : meshes())
            mesh.info();

        for (const auto& domain : domains())
            domain.info(is_outermost(domain));

        if (verbose) {
            for (const auto& vertex : vertices())
                std::cout << "[" << vertex << "] = " << vertex.index() << std::endl;

            for (const auto& mesh : meshes())
                for (const auto& triangle : mesh.triangles())
                    std::cout << "[[" << triangle.vertex(0) << "] , [" << triangle.vertex(1) << "] , ["<< triangle.vertex(2) << "]] \t = " << triangle.index() << std::endl;
        }
    }

    const Interface& Geometry::interface(const std::string& id) const {
        for (const auto& domain : domains())
            for (const auto& boundary : domain.boundaries())
                if (boundary.interface().name()==id)
                    return boundary.interface();

        // Should never append
        warning(std::string("Geometry::interface: Interface id/name \"")+id+std::string("\" not found."));
        throw OpenMEEG::BadInterface(id);
    }

    const Domain& Geometry::domain(const Vect3& p) const {
        for (const auto& domain : domains())
            if (domain.contains(p))
                return domain;

        // Should never append

        throw OpenMEEG::BadDomain("Impossible");
    }

    const Domain& Geometry::domain(const std::string& name) const {
        for (const auto& domain : domains())
            if (domain.name()==name)
                return domain;

        // Should never happen

        warning(std::string("Geometry::domain: Domain id/name \"") + name + std::string("\" not found."));
        throw OpenMEEG::BadDomain(name);
    }

    void Geometry::save(const std::string& filename) const {
        GeometryIO* io = GeometryIO::create(filename);
        try {
            io->save(*this);
        } catch (OpenMEEG::Exception& e) {
            std::cerr << e.what() << " in the file " << filename << std::endl;
            exit(e.code());
        } catch (...) {
            std::cerr << "Could not read the geometry file: " << filename << std::endl;
            exit(1);
        }
    }

    void Geometry::read_geometry_file(const std::string& filename) {
        GeometryIO* io = GeometryIO::create(filename);
        try {
            io->load(*this);
        } catch (OpenMEEG::Exception& e) {
            std::cerr << e.what() << " in the file " << filename << std::endl;
            exit(e.code());
        } catch (...) {
            std::cerr << "Could not read the geometry file: " << filename << std::endl;
            exit(1);
        }
    }

    void Geometry::import(const MeshList& mesh_list) {

        struct MeshDescription {
            std::string name;
            MeshIO*     io;
        };

        clear();

        // First read all the vertices

        std::vector<MeshDescription> mesh_descriptions;
        for (const auto& desc : mesh_list) {
            const std::string& name = desc.first;
            const std::string& path = desc.second;
            MeshIO* io = MeshIO::create(path);
            io->open();
            io->load_points(*this); 
            mesh_descriptions.push_back({ name, io });
        }

        // Second really load the meshes

        meshes().reserve(mesh_descriptions.size());
        for (const auto& desc : mesh_descriptions) {
            Mesh& mesh = add_mesh(desc.name);
            desc.io->load_triangles(mesh);
            mesh.update(true);
        }

        for (auto& desc : mesh_descriptions)
            delete desc.io;
    }

    void Geometry::read_conductivity_file(const std::string& filename) {
        try {
            typedef Utils::Properties::Named<std::string,Conductivity<double>> HeadProperties;
            HeadProperties properties(filename.c_str());

            // Store the internal conductivity of the external boundary of domain i
            // and store the external conductivity of the internal boundary of domain i

            for (auto& domain : domains()) {
                try {
                    const Conductivity<double>& cond = properties.find(domain.name());
                    domain.set_conductivity(cond.sigma());
                } catch (const Utils::Properties::UnknownProperty<HeadProperties::Id>& e) {
                    throw OpenMEEG::BadDomain(domain.name());
                }
            }
        } catch (OpenMEEG::Exception& e) {
            std::cerr << e.what() << " in the file " << filename << std::endl;
            exit(e.code());
        } catch (...) {
            std::cerr << "Could not read the conductivity file: " << filename << std::endl;
            exit(1);
        }
        conductivities = true;
    }

    // This generates unique indices for vertices and triangles which will correspond to our unknowns.

    void Geometry::generate_indices(const bool OLD_ORDERING) {

        // Either unknowns (potentials and currents) are ordered by mesh (i.e. V_1, p_1, V_2, p_2, ...) (this is the OLD_ORDERING)
        // or by type (V_1, V_2, V_3 .. p_1, p_2...) (by DEFAULT)
        // or by the user himself encoded into the vtp file.
        // if you use OLD_ORDERING make sure to iterate only once on each vertex: not to overwrite index (meshes have shared vertices).

        unsigned index = 0;
        if (!OLD_ORDERING)
            for (auto& vertex : vertices())
                vertex.index() = (invalid_vertices_.count(vertex)==0) ? index++ : unsigned(-1);

        for (auto& mesh : meshes()) {
            if (OLD_ORDERING) {
                om_error(is_nested()); // OR non nested but without shared vertices
                for (const auto& vertex : mesh.vertices())
                    vertex->index() = index++;
            }
            if (!mesh.isolated() && !mesh.current_barrier())
                for (auto& triangle : mesh.triangles())
                    triangle.index() = index++;
        }

        // even the last surface triangles (yes for EIT... )

        nb_current_barrier_triangles_ = 0;
        for (auto& mesh : meshes())
            if (mesh.current_barrier()) {
                if (!mesh.isolated()) {
                    nb_current_barrier_triangles_ += mesh.triangles().size();
                    for (auto& triangle : mesh.triangles())
                        triangle.index() = index++;
                } else {
                    for (auto& triangle : mesh.triangles())
                        triangle.index() = unsigned(-1);
                }
            }

        num_params = index;
    }

    /// \return the difference of conductivities of the 2 domains.

    double Geometry::conductivity_difference(const Mesh& m) const {
        const DomainsReference& doms = domains(m);
        double res = 0.0;
        for (const auto& domainptr : doms)
            res += domainptr->conductivity()*domainptr->mesh_orientation(m);
        return res;
    }

    /// \return 0. for non communicating meshes, 1. for same oriented meshes, -1. for different orientation

    int Geometry::oriented(const Mesh& m1,const Mesh& m2) const {
        const DomainsReference& doms = common_domains(m1,m2); // 2 meshes have either 0, 1 or 2 domains in common
        return (doms.size()==0) ? 0 : ((doms[0]->mesh_orientation(m1)==doms[0]->mesh_orientation(m2)) ? 1 : -1);
    }

    bool Geometry::selfCheck() const {

        bool OK = true;

        // Test that all meshes are well oriented and not (self-)intersecting

        for (Meshes::const_iterator mit1=meshes().begin();mit1!=meshes().end();++mit1) {
            const Mesh& mesh1 = *mit1;
            if (!mesh1.has_correct_orientation())
                warning(std::string("A mesh does not seem to be properly oriented"));

            if (mesh1.has_self_intersection()) {
                warning(std::string("Mesh is self intersecting !"));
                mesh1.info();
                OK = false;
                std::cout << "Self intersection for mesh \"" << mesh1.name() << "\"" << std:: endl;
            }
            if (is_nested()) {
                for (Meshes::const_iterator mit2=mit1+1;mit2!=meshes().end();++mit2) {
                    const Mesh& mesh2 = *mit2;
                    if (mesh1.intersection(mesh2)) {
                        warning(std::string("2 meshes are intersecting !"));
                        mesh1.info();
                        mesh2.info();
                        OK = false;
                    }
                }
            }
        }
        return OK;
    }

    bool Geometry::check(const Mesh& m) const {
        bool OK = true;
        if (m.has_self_intersection()) {
            warning(std::string("Mesh is self intersecting !"));
            m.info();
            OK = false;
        }
        for (const auto& mesh : meshes())
            if (mesh.intersection(m)) {
                warning(std::string("Mesh is intersecting with one of the mesh in geom file !"));
                mesh.info();
                OK = false;
            }
        return OK;
    }

    //  TODO: find a better name check_dipoles_in_inner_layer ?
    //  TODO: Introduce a Dipoles class.

    bool Geometry::check_inner(const Matrix& mat) const {

        if (!is_nested()) {
            std::cerr << "Dipoles are only allowed when geometry is nested." << std::endl;
            return false;
        };

        const Interface& interface = innermost_interface();
        unsigned n_outside = 0;
        for (int i=0; i<mat.nlin(); i++)
            if (!interface.contains(Vect3(mat(i,0),mat(i,1),mat(i,2))))
                ++n_outside;
        if (n_outside!=0) {
            std::cerr << n_outside << " points are outside of the inner compartment." << std::endl;
            return false;
        }
        return true;
    }

    /// Determine whether the geometry is nested or not.

    void Geometry::check_geometry_is_nested() {
        // The geometry is considered non nested if:
        // (at least) one domain is defined as being outside two or more interfaces OR....

        nested = true;
        for (const auto& domain : domains()) {
            unsigned out_interface = 0;
            if (!is_outermost(domain))
                for (const auto& boundary : domain.boundaries())
                    if (boundary.inside())
                        out_interface++;
            if (out_interface>=2) {
                nested = false;
                return;
            }
        }

        // ... if 2 interfaces are composed by a same mesh oriented into two different directions.

        for (const auto& mesh : meshes()) {
            unsigned m_oriented = 0;
            for (const auto& domain : domains())
                for (const auto& boundary : domain.boundaries())
                    for (const auto& oriented_mesh : boundary.interface().oriented_meshes())
                        if (oriented_mesh.mesh()==mesh)
                            m_oriented += oriented_mesh.orientation();
            if (m_oriented==0) {
                nested = false;
                return;
            }
        }
    }

    //  Create the vector of pairs of communicating meshes.

    void Geometry::make_mesh_pairs() {
        for (const auto& mesh1 : meshes())
            if (!mesh1.isolated())
                for (const auto& mesh2 : meshes()) {
                    const int orientation = oriented(mesh1,mesh2);
                    if ((!mesh2.isolated()) && (sigma(mesh1,mesh2)!=0.0) && orientation!=0) {
                        // Communicating meshes are used for the definition of a common domain
                        meshpairs.push_back(MeshPair(mesh1,mesh2,orientation));
                    }

                    //  Lopp only over oriented pairs of meshes.

                    if (&mesh2==&mesh1)
                        break;
                }
    }

    // Mark all meshes which touch domains with 0 conductivity

    void Geometry::mark_current_barriers() {

        // Find meshes that are current_barriers (that touch a zero-conductivity domain).
        // TODO: Should we also remove meshes for which there is no conductivity jump...
        // TODO: Instead of marking meshes and vertices, remove them from the model ?
        // TODO: isolated is not the proper name. immersed, redundant, unnecessary ?

        for (auto& domain : domains()) {
            if (almost_equal(domain.conductivity(),0.0)) {
                for (auto& boundary : domain.boundaries()) {
                    for (auto& oriented_mesh : boundary.interface().oriented_meshes()) {
                        const bool fully_immersed = (oriented_mesh.mesh().current_barrier()==true);
                        oriented_mesh.mesh().current_barrier() = true;
                        if (fully_immersed) {
                            oriented_mesh.mesh().isolated()  = true;
                            oriented_mesh.mesh().outermost() = false;
                            std::cout << "Mesh \"" << oriented_mesh.mesh().name()
                                      << "\" will be excluded from computation because it touches non-conductive domains on both sides."
                                      << std::endl;

                            //  Add all of its vertices to invalid_vertices

                            for (const auto& vertex : oriented_mesh.mesh().vertices())
                                invalid_vertices_.insert(*vertex);
                        }
                    }
                }
            }
        }

        //  Redefine outermost interface.
        //  The inside of a 0-cond domain is considered as a new outermost
        //  TODO: Can we merge this loop in the previous one ?

        for (auto& domain : domains())
            if (almost_equal(domain.conductivity(),0.0))
                for (auto& boundary : domain.boundaries())
                    if (!boundary.inside())
                        for (auto& oriented_mesh : boundary.interface().oriented_meshes())
                            if (oriented_mesh.mesh().current_barrier() && !oriented_mesh.mesh().isolated())
                                oriented_mesh.mesh().outermost() = true;

        //  Do not invalidate vertices of isolated meshes if they are shared by non isolated meshes.

        std::set<Vertex> shared_vertices;
        for (const auto& vertex : invalid_vertices_)
            for (const auto& mesh : meshes())
                if (!mesh.isolated()) {
                    const auto comp = [vertex](const Vertex* v) { return *v==vertex; };
                    const auto& vfind = std::find_if(mesh.vertices().begin(),mesh.vertices().end(),comp);
                    if (vfind!=mesh.vertices().end())
                        shared_vertices.insert(**vfind); //a shared vertex is found
                }

        for (const auto& vertex : shared_vertices)
            invalid_vertices_.erase(vertex);

        // Find the various components in the geometry.
        // The various components are separated by zero-conductivity domains.

        // Figure out the connectivity of meshes

        std::set<const Mesh*> mesh_indices;
        for (const auto& mesh : meshes())
            mesh_indices.insert(&mesh);

        std::set<const Mesh*> connected_meshes;

        while (!mesh_indices.empty()) {
            std::vector<const Mesh*> conn;
            const Mesh* se = *(mesh_indices.begin());
            mesh_indices.erase(se);
            conn.push_back(se);
            connected_meshes.insert(se);
            for (unsigned iit = 0;iit<conn.size();++iit)
                for (auto& mesh : meshes())
                    if (common_domains(*conn[iit],mesh).size()!=0 && connected_meshes.insert(&mesh).second) {
                        conn.push_back(&mesh);
                        mesh_indices.erase(&mesh);
                    }

            //  Only consider connected sets with more than one component and that are not isolated.

            if (conn.size()>1 && !conn.front()->isolated())
                independant_parts.push_back(conn);
        }

        //  Report isolated geometries

        if (independant_parts.size()>1) {
            std::cout << "The geometry is cut into several unrelated parts by non-conductive domains." << std::endl
                      << "The computation will continue. But note that electric potentials from different parts are not comparable."
                      << std::endl;

            unsigned p =0;
            for (const auto& part : independant_parts) {
                std::cout << "Part " << ++p << " is formed by meshes: { ";
                for (const auto& meshptr : part)
                    std::cout << "\"" << meshptr->name() << "\" ";
                std::cout << "}." << std::endl;
            }
        }
    }
}
