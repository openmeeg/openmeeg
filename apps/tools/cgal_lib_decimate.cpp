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

#include <cgal_lib.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

// typedef to mesh
// Domain
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domainP;
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domainP>::type TrP;
typedef CGAL::Mesh_complex_3_in_triangulation_3<TrP> C3t3P;
// Criteria
typedef CGAL::Mesh_criteria_3<TrP> Mesh_criteriaP;

// Sizing field
struct Planes
{
    typedef ::FT FT;

    typedef Mesh_domainP::Index Index;

    Planes(const OpenMEEG::Matrix &_mat, double _fs): mat(_mat), fs(_fs) {}

    FT operator()(const Point_3& p, const int, const Index&) const
    {
        bool inside = true;
        for ( unsigned i = 0; i < mat.nlin(); ++i) {
            OpenMEEG::Vect3 v(p.x()-mat(i,0), p.y() - mat(i,1), p.z() - mat(i,2));
            if ( (v(0)*mat(i,3)+v(1)*mat(i,4)+v(2)*mat(i,5) < 0) ) {
                inside = false;
            }
        }
        return (inside)?fs*1./3.:fs;
    }
    const OpenMEEG::Matrix mat;
    double fs;
};

struct Spheres
{
    typedef ::FT FT;
    typedef Mesh_domainP::Index Index;

    Spheres(const OpenMEEG::Matrix &_mat, double _fs): mat(_mat), fs(_fs) {}

    FT operator()(const Point_3& p, const int, const Index&) const
    {
        bool inside = false;
        for ( unsigned i = 0; i < mat.nlin(); ++i) {
            OpenMEEG::Vect3 v(p.x()-mat(i,0), p.y() - mat(i,1), p.z() - mat(i,2));
            if ( v.norm() < mat(i,3) ) {
                inside = true;
            }
        }
        return (inside)?fs*1./3.:fs;
    }
    const OpenMEEG::Matrix mat;
    double fs;
};

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

namespace OpenMEEG {
    class Build_Mesh2Polyhedron : public CGAL::Modifier_base<HDS> {
        const Mesh *m;
    public:
        Build_Mesh2Polyhedron(const Mesh& m): m(&m) {}
        void operator()(HDS& target) {
            CGAL::Polyhedron_incremental_builder_3<HDS> builder(target, true);
            builder.begin_surface(m->nb_vertices(), m->nb_triangles(), 6*m->nb_vertices()-12);
            std::map<const Vertex *, unsigned> map;
            unsigned i = 0;
            for ( Mesh::const_vertex_iterator vit = m->vertex_begin(); vit != m->vertex_end(); ++vit, ++i) {
                builder.add_vertex(Polyhedron::Point_3((*vit)->x(), (*vit)->y(), (*vit)->z()));
                map[*vit] = i;
            }
            for ( Mesh::const_iterator tit = m->begin(); tit != m->end(); ++tit) {
                builder.begin_facet();
                builder.add_vertex_to_facet(map[&(tit->s1())]);
                builder.add_vertex_to_facet(map[&(tit->s2())]);
                builder.add_vertex_to_facet(map[&(tit->s3())]);
                builder.end_facet();
            }
            builder.end_surface();
        }
    };

    Polyhedron Mesh2Polyhedron(const Mesh& m)
    {
        Polyhedron p;
        Build_Mesh2Polyhedron modifier(m);
        p.delegate(modifier);

        return p;
    }

    // refine a mesh
    Mesh cgal_refine(const Mesh& m_in, double radius_bound, double distance_bound, const char* sizing_field)
    {
        // Mesh criteria
        Mesh_criteriaP * criteria;
        if ( sizing_field ) {
            Matrix field(sizing_field);
            if ( field.ncol() == 6 ) {
                Planes planes(field, radius_bound);
                criteria = new Mesh_criteriaP(facet_angle=30, facet_size=planes, facet_distance=distance_bound);
            } else if ( field.ncol() == 4 ) {
                Spheres spheres(field, radius_bound);
                criteria = new Mesh_criteriaP(facet_angle=30, facet_size=spheres, facet_distance=distance_bound);
            } else {
                std::cerr << "Error: file should contain either 4 or 6 columns" << std::endl;
            }
        } else {
            criteria = new Mesh_criteriaP(facet_angle=30, facet_size=radius_bound, facet_distance=distance_bound);
        }

        // Create input polyhedron
        Polyhedron polyhedron = Mesh2Polyhedron(m_in);

        // Create domain
        Mesh_domainP domain(polyhedron);

        C3t3P c3t3;
        unsigned i = 0;
        unsigned nb_initial_points = 5;
        for ( Polyhedron::Vertex_iterator it = polyhedron.vertices_begin(); it != polyhedron.vertices_end(); it++, i++) {
            if ( i% (m_in.nb_vertices() / (1+nb_initial_points)) == 0 ) {
                c3t3.triangulation().insert(it->point());
            }
        }
        // Meshing
        CGAL::refine_mesh_3(c3t3, domain, *criteria, no_exude(), no_perturb());

        return CGAL_to_OM(c3t3);
    }

    /// decimate safely a mesh
    namespace SMS = CGAL::Surface_mesh_simplification;
    Mesh cgal_decimate(const Mesh& m_in, unsigned nb_points) {
        Polyhedron polyhedron = Mesh2Polyhedron(m_in);
        int nb_edges = 3*nb_points-6;

        std::cout << "Target: less than " << nb_points << " vertices <=> " << nb_edges << " edges.\n";

        SMS::Count_stop_predicate<Polyhedron> stop(nb_edges);

        int r = SMS::edge_collapse(
                polyhedron,
                stop,
                CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, polyhedron)).
                halfedge_index_map(get(CGAL::halfedge_external_index, polyhedron)).
                get_cost(SMS::Edge_length_cost<Polyhedron>()).
                get_placement(SMS::Midpoint_placement<Polyhedron>()));

        std::cout << "polyhedron.size_of_vertices() = " <<  polyhedron.size_of_vertices() << std::endl;
        std::cout << "polyhedron.size_of_facets() = " <<  polyhedron.size_of_facets() << std::endl;

        unsigned inum = 0;
        // write the output mesh
        Mesh m(polyhedron.size_of_vertices(), polyhedron.size_of_facets());
        std::map< Polyhedron::Vertex_const_handle, unsigned> V;
        for ( Polyhedron::Vertex_iterator vit = polyhedron.vertices_begin(); vit != polyhedron.vertices_end(); ++vit)
        {
            const Polyhedron::Vertex::Point& p = vit->point();
            m.add_vertex(Vertex(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())));
            V[vit] = inum++;
        }

        for ( Polyhedron::Facet_iterator fit = polyhedron.facets_begin(); fit != polyhedron.facets_end(); ++fit) {
            Polyhedron::Facet::Halfedge_around_facet_circulator j = fit->facet_begin();
            // Facets in polyhedral surfaces are triangles.
            CGAL_assertion( CGAL::circulator_size(j) == 3);
            const int index1 = V[j->vertex()];
            j++;
            const int index2 = V[j->vertex()];
            j++;
            const int index3 = V[j->vertex()];
            Triangle t(m.vertices()[index1], m.vertices()[index2], m.vertices()[index3]);
            m.push_back(t);
        }

        m.update();
        m.correct_global_orientation();

        std::cout << "Finished...\n" << r << " edges removed.\n" << m.nb_vertices() << " final nb_points.\n";

        return m;
    }
}
