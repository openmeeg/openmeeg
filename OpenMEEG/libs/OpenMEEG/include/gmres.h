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

#include "vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

#include <OpenMEEG_Export.h>

namespace OpenMEEG {

    // ===================================
    // = Define a Jacobi preconditionner =
    // ===================================
    template <typename M>
    class Jacobi {
    public:
        Jacobi (const M& m): J(m.nlin()) { 
            for ( unsigned i = 0; i < m.nlin(); ++i) {
                J(i, i) = 1.0 / m(i,i);
            }
        }

        Vector operator()(const Vector& g) const {
            return J*g;
        }
    
        ~Jacobi () {};
    private:
        SparseMatrix J; // diagonal
    };

    // =========================
    // = Define a GMRes solver =
    // =========================

    void GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn)
    {
        if (dy == 0.0) {
            cs = 1.0;
            sn = 0.0;
        } else if (std::abs(dy) > std::abs(dx)) {
            double temp = dx / dy;
            sn = 1.0 / sqrt( 1.0 + temp*temp );
            cs = temp * sn;
        } else {
            double temp = dy / dx;
            cs = 1.0 / sqrt( 1.0 + temp*temp );
            sn = temp * cs;
        }
    }

    void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
    {
        double temp  =  cs * dx + sn * dy;
        dy = -sn * dx + cs * dy;
        dx = temp;
    }

    template<class T>
    void Update(Vector &x, int k, T &h, Vector &s, Vector v[])
    {
            Vector y(s);
            // Backsolve:  
            for (int i = k; i >= 0; i--) {
                    y(i) /= h(i,i);
                    for (int j = i - 1; j >= 0; j--)
                            y(j) -= h(j,i) * y(i);
            }
            for (int j = 0; j <= k; j++) {
                    x += v[j] * y(j);
            }
    }

    // code taken from http://www.netlib.org/templates/cpp/gmres.h and modified
    template<class T,class P> // T should be a linear operator, and P a preconditionner
    unsigned GMRes(const T& A, const P& M, Vector &x, const Vector& b, int max_iter, double tol,unsigned m) {

        // m is the size of the Krylov subspace, if m<A.nlin(), it is a restarted GMRes (for saving memory)
        Matrix H(m+1,m);
        x.set(0.0);

        double resid;
        int i, j = 1, k;
        Vector s(m+1), cs(m+1), sn(m+1), w;

        double normb = (M(b)).norm();//(M*b).norm()
        Vector r = M(b-A*x);//M.solve(b - A * x);
        double beta = r.norm();

        if (normb == 0.0)
            normb = 1;

        if ((resid = r.norm() / normb) <= tol) {
            tol = resid;
            max_iter = 0;
            return 0;
        }
        Vector *v = new Vector[m+1];

        while (j <= max_iter) {
            v[0] = r * (1.0 / beta);
            s.set(0.0);
            s(0) = beta;

            for (i = 0; i < m && j <= max_iter; i++, j++) {
                w = M(A*v[i]); //M.solve(A * v[i]);
                for (k = 0; k <= i; k++) {
                    H(k, i) = w*v[k];
                    w -= H(k, i) * v[k];
                }
                H(i+1, i) = w.norm();
                v[i+1] = (w / H(i+1, i));

                for (k = 0; k < i; k++)
                    ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));

                GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
                ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
                ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));

                if ((resid = std::abs(s(i+1)) / normb) < tol) {
                    Update(x, i, H, s, v);
                    tol = resid;
                    max_iter = j;
                    // std::cout<<max_iter <<std::endl;
                    delete [] v;
                    return 0;
                }
            }
            Update(x, i - 1, H, s, v);
            r = M(b-A*x);//M.solve(b - A * x);
            beta = r.norm();
            if ((resid = beta / normb) < tol) {
                tol = resid;
                max_iter = j;
                // std::cout<<max_iter <<std::endl;
                delete [] v;
                return 0;
            }
        }

        tol = resid;
        delete [] v;
        return 1;
    }
}
