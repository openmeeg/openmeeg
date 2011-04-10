/* OpenMEEG

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

#if WIN32
#define _USE_MATH_DEFINES
#endif

#include <math.h>

#include "matrix.h"
#include "symmatrix.h"
#include "geometry.h"
#include "operators.h"
#include "assemble.h"

namespace OpenMEEG
{

template<class T>
void deflat(T &M, int start, int end, double coef)
{
    // deflate the Matrix
    for(int i=start; i<=end; i++) {
#ifdef USE_OMP
        #pragma omp parallel for
#endif
        for(int j=i; j<=end; j++) {
            M(i, j)+=coef;
        }
    }
}

void assemble_HM(const Geometry &geo, SymMatrix &mat, const int gauss_order)
{
    int offset = 0;

    mat = SymMatrix(geo.size());
    mat.set(0.0);

    if (geo.nb()==1) {
        operatorN(geo.getM(0), geo.getM(0), mat, 0, 0, gauss_order, offset, offset);
    }

    for (int c=0; c<geo.nb()-1; c++) {
        int offset0 = offset;
        int offset1 = offset+geo.getM(c).nbPts();
        int offset2 = offset+geo.getM(c).nbPts() + geo.getM(c).nbTrgs();
        int offset3 = offset+geo.getM(c).nbPts() + geo.getM(c).nbTrgs() + geo.getM(c+1).nbPts();

        // Computing S block first because it's needed for the corresponding N block
        if (c == 0) {
            operatorS(geo.getM(c), geo.getM(c), mat, offset1, offset1, gauss_order);
        }
        operatorS(geo.getM(c+1), geo.getM(c), mat, offset3, offset1, gauss_order);
        operatorS(geo.getM(c+1), geo.getM(c+1), mat, offset3, offset3, gauss_order);

        // Computing N block
        if (c == 0) {
            operatorN(geo.getM(c), geo.getM(c), mat, offset0, offset0, gauss_order, offset1, offset1);
        }
        operatorN(geo.getM(c+1), geo.getM(c), mat, offset2, offset0, gauss_order, offset3, offset1);
        operatorN(geo.getM(c+1), geo.getM(c+1), mat, offset2, offset2, gauss_order, offset3, offset3);

        // Computing D block
        if (c == 0) {
            operatorD(geo.getM(c), geo.getM(c), mat, offset1, offset0, gauss_order);
        }
        if (c != geo.nb()-2) {
            operatorD(geo.getM(c+1), geo.getM(c), mat, offset3, offset0, gauss_order);
        }
        operatorD(geo.getM(c), geo.getM(c+1), mat, offset1, offset2, gauss_order);
        if (c != geo.nb()-2) {
            operatorD(geo.getM(c+1), geo.getM(c+1), mat, offset3, offset2, gauss_order);
        }

        offset = offset2;
    }

    //Block multiplications
    //Because only half the Matrix is stored, only the lower part of the Matrix is treated
    offset=0;
    double K = 1.0 / (4.0 * M_PI);
    for(int c=0; c < geo.nb()-1; c++) {
        int offset0 = offset;
        int offset1 = offset0 + geo.getM(c).nbPts();
        int offset2 = offset1 + geo.getM(c).nbTrgs();
        int offset3 = offset2 + geo.getM(c+1).nbPts();
        int offset4 = offset3 + geo.getM(c+1).nbTrgs();

        //Each operator is scaled with the appropriate constant

        //Column 1
        if(c==0) {
            mult(mat, offset0, offset0, offset1, offset1, (geo.sigma_in(c)+geo.sigma_out(c))*K);
        }
        if(c==0) {
            mult(mat, offset1, offset0, offset2, offset1, -2.0*K);
        }
        mult(mat, offset2, offset0, offset3, offset1, (-geo.sigma_out(c))*K);
        mult(mat, offset3, offset0, offset4, offset1, K);

        //Column 2
        if(c==0) {
            mult(mat, offset1, offset1, offset2, offset2, K / geo.sigma_in(c) + K / geo.sigma_out(c));
        }
        mult(mat, offset2, offset1, offset3, offset2, K);
        mult(mat, offset3, offset1, offset4, offset2, -K / geo.sigma_out(c));

        //Column 3
        mult(mat, offset2, offset2, offset3, offset3, (geo.sigma_in(c+1) + geo.sigma_out(c+1))*K);
        mult(mat, offset3, offset2, offset4, offset3, -2.0*K);

        //Column 4
        mult(mat, offset3, offset3, offset4, offset4, K / geo.sigma_in(c+1) + K / geo.sigma_out(c+1));

        offset=offset2;
    }
    if (geo.nb()==1) {
        mult(mat, 0, 0, geo.getM(0).nbPts(), geo.getM(0).nbPts(), geo.sigma_in(0)*K);
    }

    // Deflate the last diagonal block of new 'mat' :
    int newsize = geo.size() - geo.getM(geo.nb() - 1).nbTrgs();
    offset = newsize - geo.getM(geo.nb() - 1).nbPts();
    deflat(mat, offset, newsize-1, mat(offset, offset) / (newsize - offset));

    mat = mat.submat(0, newsize-1);
}

void assemble_Surf2Vol(const Geometry& geo, Matrix& mat, Matrix& points)
{
    double K = 1.0 / (4.0 * M_PI);

    if (points.ncol() != 4) {
        std::cerr << " Surf2Vol: Warning no domain information provided! ";
        std::cerr << "Points are considered to be in the inner volume" << std::endl;
    }

    // Count number of points per domains

    std::vector<int> nb_pts_per_dom(geo.nb(), 0);
    for (unsigned i=0; i<points.nlin(); i++) {
        if (points.ncol() == 3) {
            ++nb_pts_per_dom[0];
        } else {
            int domain = points(i, 3);
            if (domain >= geo.nb()) {
                std::cerr << " Surf2Vol: Point " << points.getlin(i);
                std::cerr << " is outside the head. Point is considered to be in the scalp." << std::endl;
                points(i, 3) = domain = geo.nb() - 1;
            }
            ++nb_pts_per_dom[domain];
        }
    }

    std::vector<Matrix> vect_PtsInDom(geo.nb());
    for (int c=0; c<geo.nb(); ++c) {
        vect_PtsInDom[c] = Matrix(nb_pts_per_dom[c], 3);
        int iptd=0;
        for (unsigned ipt=0; ipt<points.nlin(); ipt++) { // get the points in the domain c
            if (((points.ncol()==4) && (points(ipt, 3)==c)) || ((points.ncol()==3) && (c==0))) {
                vect_PtsInDom[c](iptd, 0) = points(ipt, 0);
                vect_PtsInDom[c](iptd, 1) = points(ipt, 1);
                vect_PtsInDom[c](iptd++, 2) = points(ipt, 2);
            }
        }
    }

    mat = Matrix(points.nlin(), geo.size() - geo.getM(geo.nb()-1).nbTrgs());
    mat.set(0.0);

    int offset = 0;
    int offsetA0 = 0;
    for (int c=0; c<geo.nb()-1; c++) {
        const int offset0 = offset;
        const int offsetA = offsetA0;
        const int offsetB = offsetA + nb_pts_per_dom[c];
        const int offset1 = offset + geo.getM(c).nbPts();
        const int offset2 = offset1 + geo.getM(c).nbTrgs();
        const int offset3 = offset2 + geo.getM(c+1).nbPts();
        const int offset4 = offset3 + geo.getM(c+1).nbTrgs();

        // compute DI, i block if necessary.
        if (c==0) {
            operatorDinternal(geo.getM(c), mat, offsetA, offset0, vect_PtsInDom[c]);
            mult(mat, offsetA, offset0, offsetB, offset1,-1.0*K);
        }

        // compute DI+1, i block.
        operatorDinternal(geo.getM(c), mat, offsetB, offset0, vect_PtsInDom[c+1]);
        mult(mat, offsetB, offset0, offsetB+nb_pts_per_dom[c+1], offset1, K);

        // compute DI+1, i+1 block
        operatorDinternal(geo.getM(c+1), mat, offsetB, offset2, vect_PtsInDom[c+1]);
        mult(mat, offsetB, offset2, offsetB+nb_pts_per_dom[c+1], offset3, -K);

        // compute SI, i block if necessary.
        if (c==0) {
            operatorSinternal(geo.getM(c), mat, offsetA, offset1, vect_PtsInDom[c]);
            mult(mat, offsetA, offset1, offsetA+nb_pts_per_dom[c], offset2, K / geo.sigma_in(c));
        }

        // compute SI+1, i block
        operatorSinternal(geo.getM(c), mat, offsetB, offset1, vect_PtsInDom[c+1]);
        mult(mat, offsetA+nb_pts_per_dom[c], offset1, offsetB+nb_pts_per_dom[c+1], offset2, -K / geo.sigma_in(c+1));

        if (c<geo.nb()-2) {
            // compute SI+1, i block
            operatorSinternal(geo.getM(c+1), mat, offsetB, offset3, vect_PtsInDom[c+1]);
            mult(mat, offsetB, offset3, offsetB+nb_pts_per_dom[c+1], offset4, K / geo.sigma_in(c+1));
        }
        offset = offset2;
        offsetA0 = offsetA + nb_pts_per_dom[c];
    }
}

HeadMat::HeadMat (const Geometry& geo, const int gauss_order)
{
    assemble_HM(geo, *this, gauss_order);
}

Surf2VolMat::Surf2VolMat (const Geometry& geo, Matrix& points)
{
    assemble_Surf2Vol(geo, *this, points);
}

} // namespace OpenMEEG
