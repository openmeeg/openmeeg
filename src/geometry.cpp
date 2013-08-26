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

#include <geometry.h>
#include <geometry_reader.h>
#include <geometry_io.h>

namespace OpenMEEG {

    const Interface& Geometry::outermost_interface() const 
    {
        for ( Domains::const_iterator dit = domain_begin(); dit != domain_end(); ++dit) {
            if ( dit->outermost() ) {
                return dit->begin()->interface();
            }
        }
        // should never append
    }

    Mesh&  Geometry::mesh(const std::string& id) 
    {
        for ( iterator mit = begin() ; mit != end(); ++mit ) {
            if ( id == mit->name() ) {
                return *mit;
            }
        }
        warning(std::string("Geometry::mesh: Error mesh id/name not found: ") + id);
        // should never append
    }

    void Geometry::info() const 
    {
        if ( is_nested_ ) {
            std::cout << "This geometry is a NESTED geometry." << std::endl;
        } else {
            int shared = -1*vertices_.size();
            for (const_iterator mit = begin(); mit != end(); ++mit) {
                shared += mit->nb_vertices();
            }
            // TODO correct those are not shared vertices, but already found vertices...
            std::cout << "This geometry is a NON NESTED geometry. (There was " << shared << " demands for adding the same points)." << std::endl;
        }

        for (const_iterator mit = begin(); mit != end(); ++mit) {
            mit->info();
        }
        for ( Domains::const_iterator dit = domain_begin(); dit != domain_end(); ++dit) {
            dit->info();
        }
    }

    const Interface& Geometry::interface(const std::string& id) const 
    {
        for ( Domains::const_iterator dit = domain_begin(); dit != domain_end(); ++dit) {
            for ( Domain::const_iterator hit = dit->begin(); hit != dit->end(); ++hit) {
                if ( hit->interface().name() == id )  {
                    return hit->interface();
                }
            }
        }
        warning(std::string("Geometry::interface: Interface id/name \"") + id + std::string("\" not found."));
        // should never append
    }

    const Domain& Geometry::domain(const Vect3& p) const 
    {
        for ( Domains::const_iterator dit = domain_begin(); dit != domain_end(); ++dit) {
            if ( dit->contains_point(p) ) {
                return *dit;
            }
        }
        // should never append
    }

    const Domain& Geometry::domain(const std::string& dname) const
    {
        for ( Domains::const_iterator dit = domain_begin(); dit != domain_end(); ++dit) {
            if ( dit->name() == dname ) {
                return *dit;
            }
        }
        // should never append
        warning(std::string("Geometry::domain: Domain id/name \"") + dname + std::string("\" not found."));
    }

    void Geometry::read(const std::string& geomFileName, const std::string& condFileName) 
    {
        GeometryReader geoR(*this);

        geoR.read_geom(geomFileName);

        if ( condFileName != "" ) {
            geoR.read_cond(condFileName);
            has_cond_ = true;
        }

        // generate the indices of our unknowns
        generate_indices();

        // print info
        info();
    }

    // this generates unique indices for vertices and triangles which will correspond to our unknowns.
    void Geometry::generate_indices() 
    {
        // Either unknowns (potentials and currents) are ordered by mesh (i.e. V_1, p_1, V_2, p_2,...) 
        // or by type (V_1,V_2,V_3 .. p_1, p_2...)
        // #define CLASSIC_ORDERING // if we use classic_ordering make sure vertex do not overwrite index.. meshes have shared vertices..
        unsigned index = 0;
        #ifndef CLASSIC_ORDERING
        for ( Vertices::iterator pit = vertex_begin(); pit != vertex_end(); ++pit, index) {
            pit->index() = index++;
        }
        #endif
        for ( iterator mit = begin(); mit != end(); ++mit) {
            #ifdef CLASSIC_ORDERING
            for ( Mesh::const_vertex_iterator vit = mit->vertex_begin(); vit != mit->vertex_end(); ++vit) {
                (*vit)->index() = index++;
            }
            #endif
            for ( Mesh::iterator tit = mit->begin(); tit != mit->end(); ++tit) {
                if ( !mit->outermost() ) {
                    tit->index() = index++;
                }
            }
        } // even the last surface triangles (yes for EIT... )
        for ( iterator mit = begin(); mit != end(); ++mit) {
            for ( Mesh::iterator tit = mit->begin(); tit != mit->end(); ++tit) {
                if ( mit->outermost() ) {
                    tit->index() = index++;
                }
            }
        }
        size_ = index;
    }

    const double Geometry::sigma(const std::string& name) const 
    {
        for ( std::vector<Domain>::const_iterator dit = domain_begin(); dit != domain_end(); ++dit) {
            if ( name == dit->name() ) {
                return dit->sigma();
            }
        }
        warning(std::string("Geometry::sigma: Domain id/name \"") + name + std::string("\" not found."));
        return 0.;
    }

    const Domains Geometry::common_domains(const Mesh& m1, const Mesh& m2) const 
    {
        std::set<Domain> sdom1;
        std::set<Domain> sdom2;
        for ( Domains::const_iterator dit = domain_begin(); dit != domain_end(); ++dit) {
            if ( dit->mesh_orientation(m1) != 0 ) {
                sdom1.insert(*dit);
            }
            if ( dit->mesh_orientation(m2) != 0 ) {
                sdom2.insert(*dit);
            }
        }
        Domains doms;
        std::set_intersection(sdom1.begin(), sdom1.end(), sdom2.begin(), sdom2.end(), std::back_inserter(doms) );
        return doms;
    }

    /// \return a function (sum, difference,...) of the conductivity(ies) of the shared domain(s).
    const double Geometry::funct_on_domains(const Mesh& m1, const Mesh& m2, const Function& f) const 
    {
        Domains doms = common_domains(m1, m2);
        double ans=0.;
        for ( Domains::iterator dit = doms.begin(); dit != doms.end(); ++dit) {
            if ( f == IDENTITY ) {
                ans += dit->sigma();
            } else if ( f == INVERSE ) {
                ans += 1./dit->sigma();
            } else if ( f == INDICATOR ) {
                ans += 1.;
            } else {
                ans = 0;
            }
        }
        return ans;
    }

    /// \return the difference of conductivities of the 2 domains.
    const double  Geometry::sigma_diff(const Mesh& m) const {
        Domains doms = common_domains(m, m); // Get the 2 domains surrounding mesh m
        double  ans  = 0.;
        for ( Domains::iterator dit = doms.begin(); dit != doms.end(); ++dit) {
            ans -= dit->sigma()*dit->mesh_orientation(m);
        }
        return ans;
    }
    
    /// \return 0. for non communicating meshes, 1. for same oriented meshes, -1. for different orientation
    const double Geometry::oriented(const Mesh& m1, const Mesh& m2) const 
    {
        Domains doms = common_domains(m1, m2); // 2 meshes have either 0, 1 or 2 domains in common
        double ans = 0.;
        if ( doms.size() == 0 ) {
            return 0.;
        } else {
            return (( doms[0].mesh_orientation(m1) == doms[0].mesh_orientation(m2) ) ? 1. : -1.);
        }
    }

    const bool Geometry::selfCheck() const
    {
        bool OK = true;

        // Test that all meshes are well oriented and not (self-)intersecting
        for ( const_iterator mit1 = begin() ; mit1 != end(); ++mit1 ) {
            if ( !mit1->has_correct_orientation() ) {
                warning(std::string("A mesh does not seem to be properly oriented"));
            }
            if ( mit1->has_self_intersection() ) {
                warning(std::string("Mesh is self intersecting !"));
                mit1->info();
                OK = false;
                std::cout << "Self intersection for mesh \"" << mit1->name() << "\"" << std:: endl;
            }
            if ( is_nested_ ) {
                for ( const_iterator mit2 = mit1+1 ; mit2 != end(); ++mit2 ) {
                    if ( mit1->intersection(*mit2) ) {
                        warning(std::string("2 meshes are intersecting !"));
                        mit1->info();
                        mit2->info();
                        OK = false;
                    }
                }
            }
        }
        return OK;
    }

    const bool Geometry::check(const Mesh& m) const 
    {
        bool OK = true;
        if ( m.has_self_intersection() ) {
            warning(std::string("Mesh is self intersecting !"));
            m.info();
            OK = false;
        }
        for ( const_iterator mit = begin() ; mit != end(); ++mit ) {
            if ( mit->intersection(m) ) {
                warning(std::string("Mesh is intersecting with one of the mesh in geom file !"));
                mit->info();
                OK = false;
            }
        }
        return OK;
    }

    void Geometry::import_meshes(const Meshes& m)
    {
        meshes_.clear();
        vertices_.clear();
        unsigned n_vert_max = 0;
        unsigned iit = 0;
        std::map<const Vertex *, Vertex *> map_vertices;

        for ( Meshes::const_iterator mit = m.begin(); mit != m.end(); ++mit) {
            n_vert_max += mit->nb_vertices();
        }
        
        vertices_.reserve(n_vert_max);
        meshes_.reserve(m.size());

        for ( Meshes::const_iterator mit = m.begin(); mit != m.end(); ++mit, ++iit) {
            meshes_.push_back(Mesh(vertices_, mit->name()));
            for ( Mesh::const_vertex_iterator vit = mit->vertex_begin(); vit != mit->vertex_end(); vit++) {
                meshes_[iit].add_vertex(**vit);
                map_vertices[*vit] = *meshes_[iit].vertex_rbegin();
            }
        }

        // Copy the triangles in the geometry.
        iit = 0;
        for ( Meshes::const_iterator mit = m.begin(); mit != m.end(); ++mit, ++iit) {
            for ( Mesh::const_iterator tit = mit->begin(); tit != mit->end(); ++tit) {
                meshes_[iit].push_back(Triangle(map_vertices[(*tit)[0]], 
                                                map_vertices[(*tit)[1]], 
                                                map_vertices[(*tit)[2]]));
            }
            meshes_[iit].update();
        }
    }
#if 0
        std::set<Vertex> set_vertices;
        std::map<const Vertex, unsigned> map_vertices;

        for ( Meshes::const_iterator mit = m.begin(); mit != m.end(); ++mit, ++iit) {
            meshes_.push_back(Mesh(vertices_, mit->name()));
            for ( Mesh::const_vertex_iterator vit = mit->vertex_begin(); vit != mit->vertex_end(); vit++) {
                std::pair<std::set<Vertex>::iterator, bool> ret = set_vertices.insert(**vit);
                if ( ret.second ) {
                    map_vertices[**vit] = vertices_.size();
                    meshes_[iit].add_vertex(**vit);
                }
            }
        }

        // Copy the triangles the geometry.
        iit = 0;
        for ( Meshes::const_iterator mit = m.begin(); mit != m.end(); ++mit, ++iit) {
            for ( Mesh::const_iterator tit = mit->begin(); tit != mit->end(); ++tit) {
                meshes_[iit].push_back(Triangle(vertices_[map_vertices[tit->s1()]], 
                                                vertices_[map_vertices[tit->s2()]], 
                                                vertices_[map_vertices[tit->s3()]]));
            }
            meshes_[iit].update();
        }
    }
#endif
}
