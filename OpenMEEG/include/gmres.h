// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

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

        Jacobi(const M& m): J(m.nlin()) { 
            for (unsigned i=0; i<m.nlin(); ++i)
                J(i,i) = 1.0/m(i,i);
        }

        ~Jacobi () { }

        Vector operator()(const Vector& g) const { return J*g; }
    
    private:

        SparseMatrix J; // diagonal
    };

    // =========================
    // = Define a GMRes solver =
    // =========================

    void GeneratePlaneRotation(const double dx,const double dy,double& cs,double& sn) {
        if (dy==0.0) {
            cs = 1.0;
            sn = 0.0;
        } else if (std::abs(dy)>std::abs(dx)) {
            const double temp = dx/dy;
            sn = 1.0/sqrt(1.0+temp*temp);
            cs = temp*sn;
        } else {
            const double temp = dy/dx;
            cs = 1.0/sqrt(1.0+temp*temp);
            sn = temp*cs;
        }
    }

    void ApplyPlaneRotation(double& dx,double& dy,const double cs,const double sn) {
        const double temp = cs*dx+sn*dy;
        dy = -sn*dx+cs*dy;
        dx = temp;
    }

    template <typename T>
    void Update(Vector &x,int k,const T& h,const Vector &s,const Vector v[]) {
        Vector y(s);
        // Backsolve:  
        for (int i=k; i>=0; --i) {
            y(i) /= h(i,i);
            for (int j=i-1; j>=0; --j)
                y(j) -= h(j,i)*y(i);
        }
        for (int j=0; j<=k; ++j)
            x += v[j]*y(j);
    }

    // Code taken from http://www.netlib.org/templates/cpp/gmres.h and modified

    template <typename T,typename P> // T should be a linear operator, and P a preconditionner
    unsigned GMRes(const T& A,const P& M,Vector &x,const Vector& b,int max_iter,double tol,unsigned m) {

        // m is the size of the Krylov subspace, if m<A.nlin(), it is a restarted GMRes (for saving memory)

        Matrix H(m+1,m);
        x.set(0.0);

        Vector r = M(b-A*x); // M.solve(b-A*x);
        double beta = r.norm();

        double normb = (M(b)).norm(); // (M*b).norm()
        if (normb==0.0)
            normb = 1;

        double resid;
        if ((resid = r.norm()/normb)<=tol) {
            tol = resid;
            max_iter = 0;
            return 0;
        }

        Vector* v = new Vector[m+1];
        Vector s(m+1),cs(m+1),sn(m+1),w;

        for (int j=1; j<=max_iter; ) {
            v[0] = r*(1.0/beta);
            s.set(0.0);
            s(0) = beta;

            for (int i=0; i<m && j<=max_iter; ++i, ++j) {
                w = M(A*v[i]); // M.solve(A*v[i]);
                for (int k=0; k<=i; ++k) {
                    H(k,i) = w*v[k];
                    w -= H(k,i)*v[k];
                }
                H(i+1,i) = w.norm();
                v[i+1] = (w/H(i+1,i));

                for (int k=0; k<i; ++k)
                    ApplyPlaneRotation(H(k,i),H(k+1,i),cs(k),sn(k));

                GeneratePlaneRotation(H(i,i),H(i+1,i),cs(i),sn(i));
                ApplyPlaneRotation(H(i,i),H(i+1,i),cs(i),sn(i));
                ApplyPlaneRotation(s(i),s(i+1),cs(i),sn(i));

                if ((resid = std::abs(s(i+1))/normb)<tol) {
                    Update(x,i,H,s,v);
                    tol = resid;
                    max_iter = j;
                    delete [] v;
                    return 0;
                }
            }
            Update(x,i-1,H,s,v);
            r = M(b-A*x); // M.solve(b-A*x);
            beta = r.norm();
            if ((resid = beta/normb)<tol) {
                tol = resid;
                max_iter = j;
                delete [] v;
                return 0;
            }
        }

        tol = resid;
        delete [] v;
        return 1;
    }
}
