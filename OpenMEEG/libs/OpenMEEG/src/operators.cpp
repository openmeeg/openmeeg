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

#include <operators.h>

namespace OpenMEEG {

    void operatorDinternal(const Mesh& m, Matrix& mat, const Vertices& points, const double& coeff)
    {
        std::cout << "INTERNAL OPERATOR D..." << std::endl;
        for ( Vertices::const_iterator vit = points.begin(); vit != points.end(); ++vit)  {
            for ( Mesh::const_iterator tit = m.begin(); tit != m.end(); ++tit) {
                _operatorDinternal(*tit, *vit, mat, coeff);
            }
        }
    }

    void operatorSinternal(const Mesh& m, Matrix& mat, const Vertices& points, const double& coeff) 
    {
        std::cout << "INTERNAL OPERATOR S..." << std::endl;
        for ( Vertices::const_iterator vit = points.begin(); vit != points.end(); ++vit)  {
            for ( Mesh::const_iterator tit = m.begin(); tit != m.end(); ++tit) {
                mat(vit->index(), tit->index()) = _operatorSinternal(*tit, *vit) * coeff;
            }
        }
    }

    // General routine for applying _operatorFerguson (see this function for further comments)
    // to an entire mesh, and storing coordinates of the output in a Matrix.
    void operatorFerguson(const Vect3& x,const Mesh& m,Matrix& mat,const unsigned& offsetI,const double& coeff)
    {
        #pragma omp parallel for
        #ifndef OPENMP_3_0
        for (int i=0;i<m.vertex_size();++i) {
            const Mesh::const_vertex_iterator vit=m.vertex_begin()+i;
        #else
        for (Mesh::const_vertex_iterator vit=m.vertex_begin();vit<m.vertex_end();++vit) {
        #endif
            Vect3 v = _operatorFerguson(x, **vit, m);
            mat(offsetI + 0, (*vit)->index()) += v.x() * coeff;
            mat(offsetI + 1, (*vit)->index()) += v.y() * coeff;
            mat(offsetI + 2, (*vit)->index()) += v.z() * coeff;
        }
    }

    void operatorDipolePotDer(const Vect3& r0,const Vect3& q,const Mesh& m,Vector& rhs,const double& coeff,const unsigned gauss_order,const bool adapt_rhs) 
    {
        static analyticDipPotDer anaDPD;

        Integrator<Vect3,analyticDipPotDer>* gauss = (adapt_rhs) ? new AdaptiveIntegrator<Vect3, analyticDipPotDer>(0.001) :
                                                                   new Integrator<Vect3, analyticDipPotDer>;

        gauss->setOrder(gauss_order);
        #pragma omp parallel for private(anaDPD)
        #ifndef OPENMP_3_0
        for (int i=0;i<m.size();++i) {
            const Mesh::const_iterator tit=m.begin()+i;
        #else
        for (Mesh::const_iterator tit=m.begin();tit<m.end();++tit) {
        #endif
            anaDPD.init(*tit, q, r0);
            Vect3 v = gauss->integrate(anaDPD, *tit);
            #pragma omp critical
            {
                rhs(tit->s1().index() ) += v(0) * coeff;
                rhs(tit->s2().index() ) += v(1) * coeff;
                rhs(tit->s3().index() ) += v(2) * coeff;
            }
        }
        delete gauss;
    }

    void operatorDipolePot(const Vect3& r0, const Vect3& q, const Mesh& m, Vector& rhs, const double& coeff, const unsigned gauss_order, const bool adapt_rhs) 
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
        #ifndef OPENMP_3_0
        for (int i=0;i<m.size();++i) {
            const Mesh::const_iterator tit=m.begin()+i;
        #else
        for (Mesh::const_iterator tit=m.begin();tit<m.end();++tit) {
        #endif
            double d = gauss->integrate(anaDP, *tit);
            #pragma omp critical
            rhs(tit->index()) += d * coeff;
        }
        delete gauss;
    }

}
