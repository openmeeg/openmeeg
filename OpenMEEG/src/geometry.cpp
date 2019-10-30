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
#include <geometry_reader.h>
#include <geometry_io.h>
#include <ciso646>

namespace OpenMEEG {

    Domain&
    Geometry::outermost_domain() {
        // Search for the outermost domain and set boolean OUTERMOST on the domain in the vector domains.
        // An outermost domain is (here) defined as the only domain outside represented by only one interface.

        for (auto& domain : domains()) {
            bool outer = true;
            for (auto& halfspace : domain)
                if (halfspace.inside()) {
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
            for (const auto& halfspace : domain)
                if (!halfspace.inside()) {
                    inner = false;
                    break;
                }
            if (inner)
                return domain.front().interface();
        }

        // Should never append as this function should only be called for nested geometries.

        warning("Geometry::innerermost_interface: Error innermost interface is not defined.");
        throw OpenMEEG::BadInterface("innermost");
    }

    const Interface&
    Geometry::outermost_interface() const {
        for (const auto& domain : domains())
            if (domain.outermost())
                return domain.front().interface();

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
            domain.info();

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
            for (const auto& halfspace : domain)
                if (halfspace.interface().name()==id)
                    return halfspace.interface();

        // Should never append
        warning(std::string("Geometry::interface: Interface id/name \"")+id+std::string("\" not found."));
        throw OpenMEEG::BadInterface(id);
    }

    const Domain& Geometry::domain(const Vect3& p) const {
        for (const auto& domain : domains())
            if (domain.contains_point(p))
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

    double Geometry::sigma(const std::string& name) const {
        try {
            const Domain& dom = domain(name);
            return dom.sigma();
        } catch(...) {
            warning(std::string("Geometry::sigma: Domain id/name \"")+name+std::string("\" not found."));
            return 0.;
        }
    }

    void Geometry::read(const std::string& geomFileName, const std::string& condFileName, const bool OLD_ORDERING) {
        // clear all first
        vertices_.clear();
        meshes_.clear();
        domains_.clear();

        GeometryReader geoR(*this);
        try {
            geoR.read_geom(geomFileName);
        } catch ( OpenMEEG::Exception& e) {
            std::cerr << e.what() << " in the file " << geomFileName << std::endl;
            exit(e.code());
        } catch (...) {
            std::cerr << "Could not read the geometry file: " << geomFileName << std::endl;
            exit(1);
        }

        if (condFileName != "") {
            try {
                geoR.read_cond(condFileName);
            } catch ( OpenMEEG::Exception& e) {
                std::cerr << e.what() << " in the file " << condFileName << std::endl;
                exit(e.code());
            } catch (...) {
                std::cerr << "Could not read the conducitvity file: " << condFileName << std::endl;
                exit(1);
            }
            has_cond_ = true;

            // mark meshes that touch the 0-cond
            mark_current_barrier();
        }

        // generate the indices of our unknowns
        generate_indices(OLD_ORDERING);

        // print info
        info();
    }

    // This generates unique indices for vertices and triangles which will correspond to our unknowns.
    void Geometry::generate_indices(const bool OLD_ORDERING) {
        // Either unknowns (potentials and currents) are ordered by mesh (i.e. V_1, p_1, V_2, p_2, ...) (this is the OLD_ORDERING)
        // or by type (V_1, V_2, V_3 .. p_1, p_2...) (by DEFAULT)
        // or by the user himself encoded into the vtp file.
        // if you use OLD_ORDERING make sure to iterate only once on each vertex: not to overwrite index (meshes have shared vertices).

        if (meshes().front().triangles().front().index()==unsigned(-1)) {
            unsigned index = 0;
            if (!OLD_ORDERING)
                for (auto& vertex : vertices())
                    if (invalid_vertices_.empty() || invalid_vertices_.count(vertex)==0) {
                        vertex.index() = index++;
                    } else {
                        vertex.index() = unsigned(-1);
                    }

            for (auto& mesh : meshes()) {
                if (OLD_ORDERING) {
                    om_error(is_nested()); // OR non nested but without shared vertices
                    for (const auto& vertex : mesh.vertices())
                        vertex->index() = index++;
                }
                if (!mesh.isolated()&&!mesh.current_barrier())
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

            size_ = index;
        } else {
            std::cout << "First vertex index: " << vertices().front().index() << std::endl; // Necessary ? TODO
            size_ = vertices_.size();
            for (const auto& mesh : meshes())
                size_ += mesh.triangles().size();
        }
    }

    const Domains Geometry::common_domains(const Mesh& m1,const Mesh& m2) const {
        std::set<Domain> sdom1;
        std::set<Domain> sdom2;
        for (const auto& domain : domains()) {
            if (domain.mesh_orientation(m1)!=0)
                sdom1.insert(domain);
            if (domain.mesh_orientation(m2)!=0)
                sdom2.insert(domain);
        }
        Domains doms;
        std::set_intersection(sdom1.begin(),sdom1.end(),sdom2.begin(),sdom2.end(),std::back_inserter(doms));
        return doms;
    }

    /// \return a function (sum, difference, ...) of the conductivity(ies) of the shared domain(s).

    double Geometry::funct_on_domains(const Mesh& m1, const Mesh& m2, const Function& f) const {
        Domains doms = common_domains(m1, m2);
        double ans = 0.;
        for (Domains::iterator dit = doms.begin(); dit != doms.end(); ++dit) {
            switch (f) {
                case IDENTITY:
                    ans += dit->sigma();
                    break;
                case INVERSE:
                    ans += 1./dit->sigma();
                    break;
                case INDICATOR:
                    ans += 1.;
                    break;
                default:
                    ans = 0;
                    break;
            }
        }
        return ans;
    }

    /// \return the difference of conductivities of the 2 domains.

    double  Geometry::sigma_diff(const Mesh& m) const {
        Domains doms = common_domains(m,m); // Get the 2 domains surrounding mesh m
        double  ans  = 0.;
        for (auto& domain : doms)
            ans += domain.sigma()*domain.mesh_orientation(m);
        return ans;
    }

    /// \return 0. for non communicating meshes, 1. for same oriented meshes, -1. for different orientation

    int Geometry::oriented(const Mesh& m1,const Mesh& m2) const {
        Domains doms = common_domains(m1,m2); // 2 meshes have either 0, 1 or 2 domains in common
        return (doms.size()==0) ? 0 : ((doms[0].mesh_orientation(m1)==doms[0].mesh_orientation(m2)) ? 1 : -1);
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

    bool Geometry::check_inner(const Matrix& mat) const {
        bool OK = true;
        if (!is_nested()) {
            std::cerr << "Dipoles are only allowed when geometry is nested." << std::endl;
            OK = false;
        };
        const Interface& interface = innermost_interface();
        int n_outside = 0;
        for (int i=0; i < mat.nlin(); i++) {
            if (!interface.contains_point(Vect3(mat(i, 0), mat(i, 1), mat(i, 2)))) {
                n_outside += 1.;
            };
        };
        if (n_outside > 0) {
            std::cerr << n_outside << " points are outside of the inner most compartment." << std::endl;
            OK = false;
        }
        return OK;
    }

    void Geometry::import_meshes(const Meshes& m) {
        meshes_.clear();
        vertices_.clear();

        // Count vertices

        unsigned n_vert_max = 0;
        for (const auto& mesh : m)
            n_vert_max += mesh.vertices().size();

        // Copy vertices and triangles in the geometry.

        vertices_.reserve(n_vert_max);
        meshes_.reserve(m.size());
        std::map<const Vertex*,Vertex*> map_vertices;
        for (const auto& mesh : m) {
            meshes_.push_back(Mesh(vertices_,mesh.name()));
            Mesh& newmesh = meshes_.back();
            for (const auto& vertex : mesh.vertices()) {
                newmesh.add_vertex(*vertex);
                map_vertices[vertex] = newmesh.vertices().back();
            }
            Triangles& triangles = newmesh.triangles();
            for (const auto& triangle : mesh.triangles())
                triangles.push_back(Triangle(map_vertices[&triangle.vertex(0)], 
                                             map_vertices[&triangle.vertex(1)],
                                             map_vertices[&triangle.vertex(2)]));
            newmesh.update();
        }
    }

    // mark all meshes which touch domains with 0 conductivity
    void Geometry::mark_current_barrier() {

        // figure out the connectivity of meshes

        std::vector<int> mesh_idx;
        for (unsigned i=0;i<meshes().size();++i)
            mesh_idx.push_back(i);

        std::vector<std::vector<int>> mesh_conn;
        std::vector<int> mesh_connected, mesh_diff;
        std::set_difference(mesh_idx.begin(),mesh_idx.end(),mesh_connected.begin(),mesh_connected.end(),std::insert_iterator<std::vector<int>>(mesh_diff,mesh_diff.end()));
        while (!mesh_diff.empty()) {
            std::vector<int> conn;
            const int se = mesh_diff[0];
            conn.push_back(se);
            mesh_connected.push_back(se);
            for (unsigned iit = 0;iit<conn.size();++iit) {
                const Mesh& me = meshes()[conn[iit]];
                for (auto& mesh : meshes()) {
                    if (not almost_equal(sigma(me,mesh),0.0)) {
                        const int id = &mesh-&meshes().front();
                        std::vector<int>::iterator ifind = std::find(mesh_connected.begin(),mesh_connected.end(),id);
                        if (ifind==mesh_connected.end()) {
                            mesh_connected.push_back(id);
                            conn.push_back(id);
                        }
                    }
                }
            }
            mesh_conn.push_back(conn);
            std::sort(mesh_connected.begin(), mesh_connected.end());
            mesh_diff.clear();
            std::set_difference(mesh_idx.begin(), mesh_idx.end(), mesh_connected.begin(), mesh_connected.end(), std::insert_iterator<std::vector<int> >(mesh_diff, mesh_diff.end()));
        }

        // find isolated meshes and touch 0-cond meshes;
        std::set<std::string> touch_0_mesh;
        for (Domains::iterator dit = domains_.begin(); dit != domains_.end(); ++dit) {
            if (almost_equal(dit->sigma(),0.0)) {
                for (Domain::iterator hit = dit->begin(); hit != dit->end(); ++hit) {
                    for (Interface::iterator omit = hit->first.begin(); omit != hit->first.end(); ++omit) {
                        omit->mesh().current_barrier() = true;
                        std::pair<std::set<std::string>::iterator, bool> ret = touch_0_mesh.insert(omit->mesh().name());
                        if (!ret.second) {
                            omit->mesh().isolated() = true;
                            omit->mesh().outermost() = false;
                            std::cout<<"Mesh \""<<omit->mesh().name()<<"\" will be excluded from computation because it touches non-conductive domains on both sides."<<std::endl;
                            //add all of its vertices to invalid_vertices
                            for (const auto& vertex : omit->mesh().vertices())
                                invalid_vertices_.insert(*vertex);
                        }
                    }
                }
            }
        }

        std::vector<std::vector<int> > new_conn;
        for (std::vector<std::vector<int> >::const_iterator git = mesh_conn.begin(); git != mesh_conn.end(); ++git) {
            if (git->size()>1||!meshes()[*git->begin()].isolated()) {
                new_conn.push_back(*git);
            }
        }
        mesh_conn = new_conn;

        //do not delete shared vertices
        std::set<Vertex> shared_vtx;
        for (std::set<Vertex>::const_iterator vit = invalid_vertices_.begin(); vit != invalid_vertices_.end(); ++vit)
            for (const auto& mesh : meshes())
                if (!mesh.isolated()) {
                    const auto comp = [vit](const Vertex* v) { return *v==*vit; };
                    const std::vector<Vertex*>::const_iterator vfind = std::find_if(mesh.vertices().begin(),mesh.vertices().end(),comp);
                    if (vfind!=mesh.vertices().end())
                        shared_vtx.insert(**vfind); //a shared vertex is found
                }

        for (std::set<Vertex>::const_iterator vit = shared_vtx.begin(); vit != shared_vtx.end(); ++vit)
            invalid_vertices_.erase(*vit);

        //redefine outermost interface
        //the inside of a 0-cond domain is considered as a new outermost

        for (auto& domain : domains())
            if (almost_equal(domain.sigma(),0.0))
                for (auto& halfspace : domain)
                    if (!halfspace.inside())
                        for (auto& interface : halfspace.first)
                            if (interface.mesh().current_barrier() && !interface.mesh().isolated())
                                interface.mesh().outermost() = true;

        //detect isolated geometries
        if (mesh_conn.size()>1) {
            std::cout<<"The geometry is cut into several unrelated parts by non-conductive domains."<<std::endl;
            std::cout<<"The computation will continue. But please note that the electric potentials from different parts are no longer comparable."<<std::endl;

            for (unsigned iit = 0, p = 0; iit < mesh_conn.size(); ++iit) {
                std::cout<<"Part "<<++p<<" is formed by meshes: { ";
                for (unsigned miit = 0; miit < mesh_conn[iit].size(); ++miit) {
                    std::cout<<"\""<<meshes()[mesh_conn[iit][miit]].name()<<"\" ";
                }
                std::cout<<"}."<<std::endl;
            }
        }
        //count geo_group
        for (unsigned git = 0; git < mesh_conn.size(); ++git) {
            Strings gg;
            for (unsigned mit = 0; mit < mesh_conn[git].size(); ++mit) {
                gg.push_back(meshes()[mesh_conn[git][mit]].name());
            }
            geo_group_.push_back(gg);
        }

    }
}
