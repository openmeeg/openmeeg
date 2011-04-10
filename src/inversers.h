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

#ifndef OPENMEEG_INVERSE_H
#define OPENMEEG_INVERSE_H

#include "matrix.h"
#include "symmatrix.h"
#include "vector.h"
#include "sparse_matrix.h"
#include "fast_sparse_matrix.h"

#include "DLLDefinesOpenMEEG.h"

#define EPSILON 1e-6
#define MINRES_TOL 1e-5

namespace OpenMEEG {

    // ========================================================
    // = Define Hessian matrices for linear inversion methods =
    // ========================================================

    class OPENMEEG_EXPORT MN_Hessian {

        const Matrix &Transfer;
        const double alpha;

    public:

        MN_Hessian(const Matrix &TransferMat, const double &Alpha):Transfer(TransferMat),alpha(Alpha) {}

        inline Vector operator * (const Vector &x) const { return Transfer.tmult(Transfer*x)+alpha*x; }
    };

    class OPENMEEG_EXPORT WMN_Hessian {

        const Matrix &Transfer;
        const double alpha;
        Vector weights;

    public:

        WMN_Hessian(const Matrix &TransferMat, const double &Alpha):Transfer(TransferMat),alpha(Alpha) {
            Vector v(Transfer.ncol());
            for(size_t i = 0; i < weights.size(); ++i) {
                Vector col = Transfer.getcol(i);
                v(i) = pow(col.norm(),2);
            }
            weights = v;
        }

        inline Vector operator * (const Vector &x) const {
            return Transfer.tmult(Transfer*x)+alpha*(weights.kmult(x));
        }

    };

    class OPENMEEG_EXPORT HEAT_Hessian {

        const Matrix &m_transfer;
        const FastSparseMatrix &m_mat;
        const FastSparseMatrix &m_mat_t;
        const double m_alpha;

    public:

        HEAT_Hessian(const Matrix &transfer,
                         const FastSparseMatrix &mat,
                         const FastSparseMatrix &mat_t,
                         const double &alpha):
        m_transfer(transfer),m_mat(mat),m_mat_t(mat_t),m_alpha(alpha) {}

        inline Vector operator * ( const Vector &x) const {
            return m_transfer.tmult(m_transfer*x)+m_alpha*(m_mat_t*(m_mat*x));
        }
    };

    // ========================================================
    template<class T> // T should be a linear operator
    size_t MinRes2(const T& A,const Vector& b,Vector& x0,double tol) {

        size_t n_max = 10000;
        size_t n = 1;
        size_t N = x0.size();
        Vector v(N);
        v.set(0.0);
        Vector v_hat = b-A*x0;
        double beta = v_hat.norm();
        Vector v_old(v.size());
        Vector Av(v.size());
        double c = 1;
        double c_old = 1;
        double s_old = 0;
        double s = 0;
        Vector w(N);
        w.set(0.0);
        Vector w_oold(N);
        Vector w_old(w,DEEP_COPY);
        double eta = beta;
        Vector xMR = x0;
        double norm_rMR = beta;
        double norm_r0 = beta;
        double c_oold,s_oold,r1_hat,r1,r2,r3,alpha,beta_old;
        while ((n < n_max+1) && (norm_rMR/norm_r0 > tol) ) {
            n = n+1;
            //Lanczos
            v_old = v;
            v = v_hat*(1.0/beta);
            Av = A*v;
            alpha = v*Av;
            v_hat = Av-alpha*v-beta*v_old;
            beta_old = beta;
            beta = v_hat.norm();
            //QR factorization
            c_oold = c_old;
            c_old = c;
            s_oold = s_old;
            s_old = s;
            r1_hat = c_old*alpha-c_oold*s_old*beta_old;
            r1  =  sqrt(r1_hat*r1_hat+beta*beta);
            r2  =  s_old*alpha+c_oold*c_old*beta_old;
            r3  =  s_oold*beta_old;
            //Givens rotation
            c = r1_hat/r1;
            s = beta/r1;
            //update
            w_oold = w_old;
            w_old = w;
            w = (v-r3*w_oold-r2*w_old)*(1.0/r1);
            xMR += c*eta*w;
            norm_rMR *= fabs(s);
            eta *= -s;
        }
        std::cout<<"\r";
        return n;
    }

    // ========================================================

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
            for (int j = 0; j <= k; j++)
                x += v[j] * y(j);
        }

    // code taken from http://www.netlib.org/templates/cpp/gmres.h and modified
    template<class T,class P> // T should be a linear operator, and P a preconditionner
    size_t GMRes(const T& A, const P& M, Vector &x, const Vector& b, int max_iter, double tol,unsigned m) {

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

    // ===========================================
    // = Define all the linear inversion methods =
    // ===========================================

    template<class T>
    void LIN_inverse (Matrix& EstimatedData, const T& hess, const Matrix& GainMatrix, const Matrix& Data) {
        size_t nT = Data.ncol();
        EstimatedData = Matrix(GainMatrix.ncol(),nT);

        #ifdef USE_OMP
        #pragma omp parallel for
        #endif
        for(long frame=0;frame<(long)nT;frame++) { // loop over frame
            Vector m = Data.getcol((size_t) frame);
            Vector v(GainMatrix.ncol()); v.set(0.0);

            //==========  Invert =======================//
            size_t niter = MinRes2(hess,GainMatrix.tmult(m),v,MINRES_TOL);

            for(size_t i=0;i<EstimatedData.nlin();i++) EstimatedData(i,(size_t) frame) = v(i);

            #ifdef USE_OMP
            #pragma omp critical
            #endif
            std::cout << ">> Frame " << frame+1 << " / " << nT
                      << " : Rel. Err. = " << (GainMatrix*v-m).norm()/m.norm()
                      << " : Nb. iter. MinRes = " << niter
                      << std::endl;
        }
    }

    void compute_mn (Matrix& EstimatedData, const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight) {
        Matrix eye(GainMatrix.nlin(),GainMatrix.nlin());
        eye.set(0);
        for (size_t i=0;i<GainMatrix.nlin();++i)
            eye(i,i) = SmoothWeight;
        EstimatedData = GainMatrix.transpose() * (GainMatrix * GainMatrix.transpose() + eye).inverse() * Data;
    }

    // ================= Iterative Mininum norm inversion =======================//

    class IMN_inverse: public Matrix {
    public:
        IMN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight);
        virtual ~IMN_inverse () {};
    };

    IMN_inverse::IMN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight) {
        std::cout << "Running Iterative MN inversion" << std::endl;
        MN_Hessian hess(GainMatrix,SmoothWeight);
        LIN_inverse(*this,hess,GainMatrix,Data);
    }

    // ================= Mininum norm inversion =======================//

    class MN_inverse: public Matrix {
    public:
        MN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight);
        virtual ~MN_inverse () {};
    };

    MN_inverse::MN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight) {
        std::cout << "Running MN inversion" << std::endl;
        compute_mn(*this,Data,GainMatrix,SmoothWeight);
    }

    // ================= Weighted Mininum norm inversion =======================//

    class WMN_inverse: public Matrix {
    public:
        WMN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight);
        virtual ~WMN_inverse () {};
    };

    WMN_inverse::WMN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight) {
        std::cout << "Running WMN inversion" << std::endl;
        WMN_Hessian hess(GainMatrix,SmoothWeight);
        LIN_inverse(*this,hess,GainMatrix,Data);
    }

    // ================= Gradient based Mininum norm inversion ================ //

    class HEAT_inverse: public Matrix {
    public:
        HEAT_inverse (const Matrix& Data, const Matrix& GainMatrix, const SparseMatrix& SmoothMatrix, double SmoothWeight);
        virtual ~HEAT_inverse () {};
    };

    HEAT_inverse::HEAT_inverse (const Matrix& Data, const Matrix& GainMatrix, const SparseMatrix& SmoothMatrix, double SmoothWeight) {
        std::cout << "Running HEAT inversion" << std::endl;
        FastSparseMatrix fastSmoothMatrix(SmoothMatrix);
        FastSparseMatrix fastSmoothMatrix_t(SmoothMatrix.transpose());
        HEAT_Hessian hess(GainMatrix,fastSmoothMatrix,fastSmoothMatrix_t,SmoothWeight);
        LIN_inverse(*this,hess,GainMatrix,Data);
    }

}

#endif /* OPENMEEG_INVERSE_H */

