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

#include <operators.h>

namespace OpenMEEG {

    void operatorDinternal(const Mesh& m, Matrix& mat, const unsigned offsetI, const unsigned offsetJ, const Matrix& points) 
    {
        std::cout<<"INTERNAL OPERATOR D..."<<std::endl;
        for ( unsigned i = offsetI; i < offsetI + points.nlin(); ++i)  {
            Vect3 pt(points(i-offsetI, 0), points(i-offsetI, 1), points(i-offsetI, 2));
            for ( Mesh::const_iterator tit = m.begin(); tit != m.end(); ++tit) {
                _operatorDinternal(*tit, pt, mat);
            }
        }
    }

    void operatorSinternal(const Mesh& m, Matrix& mat, const unsigned offsetI, const unsigned offsetJ, const Matrix& points) 
    {
        std::cout<<"INTERNAL OPERATOR S..."<<std::endl;
        for ( unsigned i=offsetI; i < offsetI + points.nlin(); ++i) {
            Vect3 pt(points(i-offsetI, 0), points(i-offsetI, 1), points(i-offsetI, 2));
            // for (Mesh::const_iterator tit = m.begin(); tit != m.end(); ++tit) {
            for ( unsigned j=offsetJ; j<offsetJ+m.nb_triangles(); ++j) {
                mat(i, j) = _operatorSinternal(m[j-offsetJ], pt);
            }
        }
    }

    // General routine for applying _operatorFerguson (see this function for further comments)
    // to an entire mesh, and storing coordinates of the output in a Matrix.
    void operatorFerguson(const Vect3& x, const Mesh& m, SparseMatrix& mat, const unsigned& offsetI, const double& coef)
    {
        #pragma omp parallel for
        for ( Mesh::const_vertex_iterator vit = m.vertex_begin(); vit < m.vertex_end(); ++vit) {
            Vect3 v = _operatorFerguson(x, **vit, m);
            mat(offsetI + 0, (*vit)->index()) += v.x() * coef;
            mat(offsetI + 1, (*vit)->index()) += v.y() * coef;
            mat(offsetI + 2, (*vit)->index()) += v.z() * coef;
        }
    }

    void operatorDipolePotDer(const Vect3& r0, const Vect3& q, const Mesh& m, Vector& rhs, const unsigned gauss_order, const bool adapt_rhs) 
    {
        static analyticDipPotDer anaDPD;

        Integrator<Vect3, analyticDipPotDer>* gauss;
        if ( adapt_rhs ) {
            gauss = new AdaptiveIntegrator<Vect3, analyticDipPotDer>(0.001);
        } else {
            gauss = new Integrator<Vect3, analyticDipPotDer>;
        }

        gauss->setOrder(gauss_order);
        #pragma omp parallel for private(anaDPD)
        for ( Mesh::const_iterator tit = m.begin(); tit < m.end(); ++tit) {
            anaDPD.init(*tit, q, r0);
            Vect3 v = gauss->integrate(anaDPD, *tit);
            #pragma omp critical
            {
                rhs(tit->s1().index() ) += v(0);
                rhs(tit->s2().index() ) += v(1);
                rhs(tit->s3().index() ) += v(2);
            }
        }
        delete gauss;
    }

    void operatorDipolePot(const Vect3& r0, const Vect3& q, const Mesh& m, Vector& rhs, const unsigned gauss_order, const bool adapt_rhs) 
    {
        static analyticDipPot anaDP;

        anaDP.init(q, r0);
        Integrator<double, analyticDipPot> *gauss;
        if ( adapt_rhs ) {
            gauss = new AdaptiveIntegrator<double, analyticDipPot>(0.001);
        } else {
            gauss = new Integrator<double, analyticDipPot>;
        }

        gauss->setOrder(gauss_order);
        #pragma omp parallel for
        for ( Mesh::const_iterator tit = m.begin(); tit < m.end(); ++tit) {
            double d = gauss->integrate(anaDP, *tit);
            #pragma omp critical
            rhs(tit->index()) += d;
        }
        delete gauss;
    }

} // namespace OpenMEEG
