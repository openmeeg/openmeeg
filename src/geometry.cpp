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

#include "geometry.h"
#include "reader.h"

namespace OpenMEEG {

    const Interface& Geometry::outermost_interface() const 
    {
        for (Domains::const_iterator dit = domain_begin(); dit != domain_end(); dit++) {
            if ( dit->outermost()) {
                return dit->begin()->interface();
            }
        }
    }

    const unsigned Geometry::nb_trianglesoutermost() const 
    {
        unsigned nb_t = 0;
        for (const_iterator mit = this->begin(); mit != this->end(); mit++) {
            if ( mit->outermost()) {
                nb_t += mit->nb_triangles();
            }
        }
        return nb_t;
    }

    const Mesh&  Geometry::mesh(const std::string &id) const 
    {
        for (const_iterator mit = this->begin() ; mit != this->end(); mit++ ) {
            if ( id == mit->name()) {
                return *mit;
            }
        }
        warning(std::string("Error mesh id/name not found: ") + id);
    }

    Mesh&  Geometry::mesh(const std::string &id) 
    {
        for (iterator mit = this->begin() ; mit != this->end(); mit++ ) {
            if ( id == mit->name()) {
                return *mit;
            }
        }
        warning(std::string("Error mesh id/name not found: ") + id);
    }

    void Geometry::info() const 
    {
        // for (const_iterator mit = this->begin(); mit != this->end(); mit++) {
            // mit->info();
        // }
        for (Domains::const_iterator dit = this->domain_begin(); dit != this->domain_end(); dit++) {
            dit->info();
        }
    }

    const Interface& Geometry::interface(const std::string& id) const 
    {
        for (Domains::const_iterator dit = this->domain_begin(); dit != this->domain_end(); dit++) {
            for (Domain::const_iterator hit = dit->begin(); hit != dit->end(); hit++) {
                if ( hit->interface().name() == id )  {
                    return hit->interface();
                }
            }
        }
        // should never append
    }

    const Domain& Geometry::domain(const Vect3& p) const 
    {
        for ( Domains::const_iterator dit = this->domain_begin(); dit != this->domain_end(); dit++) {
            if ( dit->contains_point(p) ) {
                return *dit;
            }
        }
        // should never append
    }

    const Domain& Geometry::domain(const std::string& dname) const
    {
        for (Domains::const_iterator dit = this->domain_begin(); dit != this->domain_end(); dit++) {
            if ( dit->name() == dname ) {
                return *dit;
            }
        }
        // should never append
        warning(std::string("Geometry::domain: Domain id/name \"") + dname + std::string("\" not found."));
    }

    void Geometry::read(const std::string geomFileName, const std::string condFileName) 
    {
        has_cond() = false; // default parameter

        read_geom(geomFileName);

        // generate the indices of our unknowns
        geom_generate_indices();

        if ( condFileName != "" ) {
            read_cond(condFileName);
            has_cond() = true;
        }
        has_cond() = true;

        info();
    }

    // this generates unique indices for vertices and triangles which will correspond to our unknowns.
    void Geometry::geom_generate_indices() 
    {
        // Either unknowns (potentials and currents) are ordered by mesh (i.e. V_1, p_1, V_2, p_2,...) 
        // or by type (V_1,V_2,V_3 .. p_1, p_2...)
        #define CLASSIC_ORDERING
        unsigned index = 0;
        #ifndef CLASSIC_ORDERING
        for ( Vertices::iterator pit = this->vertex_begin(); pit != this->vertex_end(); pit++, index) {
            pit->index() = index++;
        }
        #endif
        for ( iterator mit = this->begin(); mit != this->end(); mit++) {
            #ifdef CLASSIC_ORDERING
            for ( Mesh::const_vertex_iterator vit = mit->vertex_begin(); vit != mit->vertex_end(); vit++) {
                (*vit)->index() = index++;
            }
            #endif
            for ( Mesh::iterator tit = mit->begin(); tit != mit->end(); tit++) {
                tit->index() = index++;
            }
        }
        this->size() = index;
    }

    const double Geometry::sigma(const std::string& name) const 
    {
        for ( std::vector<Domain>::const_iterator d = domain_begin(); d != domain_end(); d++) {
            if ( name == d->name() ) {
                return d->sigma();
            }
        }
        return 0.; // TODO throw error unknownDomain
    }

    const Domains Geometry::common_domains(const Mesh& m1, const Mesh& m2) const 
    {
        std::set<Domain> sdom1;
        std::set<Domain> sdom2;
        for ( Domains::const_iterator dit = domain_begin(); dit != domain_end(); dit++) {
            if ( dit->meshOrient(m1) != 0 ) {
                sdom1.insert(*dit);
            }
            if ( dit->meshOrient(m2) != 0 ) {
                sdom2.insert(*dit);
            }
        }
        Domains doms;
        std::set_intersection(sdom1.begin(), sdom1.end(), sdom2.begin(), sdom2.end(), std::back_inserter(doms) );
        return doms;
    }

    // return the (sum) conductivity(ies) of the shared domain(s).
    const double Geometry::funct_on_domains(const Mesh& m1, const Mesh& m2, const char& f) const 
    {
        Domains doms = common_domains(m1, m2);
        double ans=0.;
        for ( Domains::iterator dit = doms.begin(); dit != doms.end(); dit++) {
            if ( f == '+' ) {
                ans += dit->sigma();
            } else if ( f == '-' ) {
                ans -= -1.*dit->sigma();
            } else if ( f == '/' ) {
                ans += 1./dit->sigma();
            } else {
                ans += 1.;
            }
        }
        return ans;
    }

    const double Geometry::oriented(const Mesh& m1, const Mesh& m2) const 
    {
        Domains doms = common_domains(m1, m2);
        double ans = 0.;
        if ( doms.size() == 2 ) { // TODO Maureen comment on the cylinder
            return 1.;
        } else if ( doms.size() == 1 ) {
            return (( doms[0].meshOrient(m1) == doms[0].meshOrient(m2) ) ? 1. : -1.);
        } else {
            return 0.;
        }
    }

    bool Geometry::selfCheck() const  // TODO: general enough ?
    {
        bool OK = true;
        for (const_iterator mit1 = this->begin() ; mit1 != this->end(); mit1++ ) {
            if ( !mit1->has_correct_orientation() ) {
                warning(std::string("A mesh does not seem to be properly oriented"));
            }
            if ( mit1->has_self_intersection() ) {
                warning(std::string("Mesh is self intersecting !"));
                mit1->info();
                OK = false;
                std::cout << "Self intersection for mesh \"" << mit1->name() << "\"" << std:: endl;
            }
            for ( const_iterator mit2 = mit1+1 ; mit2 != this->end(); mit2++ ) {
                if ( mit1->intersection(*mit2) ) {
                    warning(std::string("2 meshes are intersecting !"));
                    mit1->info();
                    mit2->info();
                    OK = false;
                }
            }
        }
        return OK;
    }

    bool Geometry::check(const Mesh& m) const 
    {
        bool OK = true;
        if ( m.has_self_intersection() ) {
            warning(std::string("Mesh is self intersecting !"));
            m.info();
            OK = false;
        }
        for ( const_iterator mit = this->begin() ; mit != this->end(); mit++ ) {
            if ( mit->intersection(m) ) {
                warning(std::string("Mesh is intersecting with one of the mesh in geom file !"));
                mit->info();
                OK = false;
            }
        }
        return OK;
    }
}
