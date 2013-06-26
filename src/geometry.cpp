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

    const size_t Geometry::nb_trianglesoutermost_() const {
        size_t nb_t = 0;
        for (const_iterator mit = this->begin(); mit != this->end(); mit++) {
            if (mit->outermost()) {
                nb_t += mit->nb_triangles();
            }
        }
        return nb_t;
    }

    const Mesh&  Geometry::mesh(const std::string &id) const {
        for (Meshes::const_iterator mit = this->begin() ; mit != this->end(); mit++ ) {
            if (id == mit->name()) {
                return *mit;
            }
        }
        std::cerr << "Error mesh id/name not found: " << id << std::endl;
    }

    Mesh&  Geometry::mesh(const std::string &id) {
        for (Meshes::iterator mit = this->begin() ; mit != this->end(); mit++ ) {
            if (id == mit->name()) {
                return *mit;
            }
        }
        std::cerr << "Error mesh id/name not found: " << id << std::endl;
    }

    void Geometry::info() const {
        for (Domains::const_iterator dit = this->domain_begin(); dit != this->domain_end(); dit++) {
            dit->info();
        }
    }

    const Domain Geometry::get_domain(const Vect3& p) const {

        for (Domains::const_iterator dit = this->domain_begin(); dit != this->domain_end(); dit++) {
            if (dit->contains_point(p)) {
                return *dit;
            }
        }
        return Domain(); // should never append
    }

    void Geometry::read(const char* geomFileName, const char* condFileName) {

        destroy();

        // vertices_.reserve(3000);
        has_cond() = false; // default parameter

        read_geom(geomFileName);

        // updates
        geom_generate_indices();

        for (Vertices::const_iterator vit1 = this->vertex_begin(); vit1 != this->vertex_end(); vit1++) {
            std::cout << "index1: " << vit1->index() << std::endl;
        }

        for (Geometry::const_iterator mit1 = this->begin(); mit1 != this->end(); mit1++) {
            for(Mesh::const_iterator tit = mit1->begin(); tit != mit1->end(); tit++) {
                std::cout << "indexTri: " << (tit)->index() << std::endl;
            }
            for(Vertices::const_iterator vit = mit1->all_vertices().begin(); vit != mit1->all_vertices().end(); vit++) {
                std::cout << "index2: " << (vit)->index() << std::endl;
            }
            for(Mesh::const_vertex_iterator vit = mit1->vertex_begin(); vit != mit1->vertex_end(); vit++) {
                std::cout << "index3: " << (*vit)->index() << std::endl;
            }
            std::cout << "\t\t\t\tDeuxieme" << std::endl;
        }

        if(condFileName) {
            read_cond(condFileName);
            has_cond() = true;
        }
        has_cond() = true;

        info();
    }

    void Geometry::geom_generate_indices() { // this generates unique indices for vertices and triangles which will correspond to our unknowns.

        size_t index = 0;
        for ( Vertices::iterator pit = this->vertex_begin(); pit != this->vertex_end(); pit++, index++) {
            pit->index() = index;
        }

        for (iterator mit = this->begin(); mit != this->end(); mit++) {
            for (Mesh::iterator tit = mit->begin(); tit != mit->end(); tit++, index++) {
                tit->index() = index;
            }
            for(Mesh::const_vertex_iterator vit = mit->vertex_begin(); vit != mit->vertex_end(); vit++) {
                std::cout << "index3: " << (*vit)->index() << std::endl;
            }
        }
        this->size() = index;
    }

    double Geometry::sigma(const std::string& name) const {
        for (std::vector<Domain>::const_iterator d = domain_begin(); d != domain_end(); d++) {
            if (name == d->name()) {
                return (d->sigma());
            }
        }
        return -1.; // TODO throw error unknownDomain
    }

    double Geometry::oriented(const Mesh& m1, const Mesh& m2) const {
        Domains doms = common_domains(m1, m2);
        double ans = 0.;
        if (doms.size() == 2) { // TODO Maureen comment on the cylinder
            return 1.;
        }
        return (doms[0].meshOrient(m1) == doms[0].meshOrient(m2))?1.:-1.;
    }

    bool Geometry::selfCheck() const { // TODO: something else
        return true;
        // bool OK = true;
        // for(int i = 0; i < nb(); ++i)
        // {
            // const Mesh& m1 = getM(i);
            // if (!m1.has_correct_orientation()) {
                // warning(std::string("A mesh does not seem to be properly oriented"));
            // }
            // if(m1.has_self_intersection())
            // {
                // warning(std::string("Mesh is self intersecting !"));
                // m1.info();
                // OK = false;
                // std::cout << "Self intersection for mesh number " << i << std:: endl;
            // }
            // for(int j = i+1; j < nb(); ++j)
            // {
                // const Mesh& m2 = getM(j);
                // if(m1.intersection(m2))
                // {
                    // warning(std::string("2 meshes are intersecting !"));
                    // m1.info();
                    // m2.info();
                    // OK = false;
                // }
            // }
        // }
        // return OK;
    }

    /*
    // Mesh& Geometry::interior_mesh(const std::string& name) const {
        // for (std::vector<Domain>::const_iterator d = domain_begin(); d != domain_end(); d++) {
            // if ( name == d->name() ) {
                // return (this->begin()); // TODO Correct that now
            // }
        // }
        // return -1.; // TODO throw error unknownDomain
    // }


       bool Geometry::check(const Mesh& m) const {
       bool OK = true;
       if(m.has_self_intersection())
       {
       warning(std::string("Mesh is self intersecting !"));
       m.info();
       OK = false;
       }
       for(int i = 0; i < nb(); ++i)
       {
       const Mesh& m1 = getM(i);
       if(m1.intersection(m))
       {
       warning(std::string("Mesh is intersecting with one of the mesh in geom file !"));
       m1.info();
       OK = false;
       }
       }
       return OK;
       }

       int Geometry::getDomain(const Vect3& p) const {
       for (int i=0;i<nb();++i)
       if ((this->getM(i)).contains_vertex(p))
       return i;
       return nb();
       }
       */

}
