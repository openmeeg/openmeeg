/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
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

#include "sparse_matrice_dcl.h"

vecteur sparse_matrice::operator*(const vecteur &x) const
{
    vecteur ret(nlin());
    ret.set(0);

    Tank::const_iterator it;
    for(it = m_tank.begin(); it != m_tank.end(); ++it) {
        size_t i = it->first.first;
        size_t j = it->first.second;
        double val = it->second;
        ret(i) += val * x(j);
    }

    return ret;
}

matrice sparse_matrice::operator*(const matrice &mat) const
{
    assert(ncol()==mat.nlin());
    matrice out(nlin(),mat.ncol());
    out.set(0.0);

    Tank::const_iterator it;
    for(it = m_tank.begin(); it != m_tank.end(); ++it) {
        size_t i = it->first.first;
        size_t j = it->first.second;
        double val = it->second;
        for(size_t k = 0; k < mat.ncol(); ++k) {
            out(i,k) += val * mat(j,k);
        }
    }

    return out;
}

sparse_matrice sparse_matrice::transpose() const {
    sparse_matrice tsp;
    const_iterator it;
    for(it = m_tank.begin(); it != m_tank.end(); ++it) {
        size_t i = it->first.first;
        size_t j = it->first.second;
        tsp(j,i) = it->second;
    }
    return tsp;
}

void sparse_matrice::info() const {
    if ((nlin() == 0) || (ncol() == 0) || m_tank.empty()) {
        std::cout << "Matrix Empty" << std::endl;
        return;
    }

    std::cout << "Dimensions : " << nlin() << " x " << ncol() << std::endl;

    double minv = m_tank.begin()->second;
    double maxv = m_tank.begin()->second;
    size_t mini = 0;
    size_t maxi = 0;
    size_t minj = 0;
    size_t maxj = 0;

    Tank::const_iterator it;
    for(it = m_tank.begin(); it != m_tank.end(); ++it) {
            if (minv > it->second) {
                minv = it->second;
                mini = it->first.first;
                minj = it->first.second;
            } else if (maxv < it->second) {
                maxv = it->second;
                maxi = it->first.first;
                maxj = it->first.second;
            }
    }

    std::cout << "Min Value : " << minv << " (" << mini << "," << minj << ")" << std::endl;
    std::cout << "Max Value : " << maxv << " (" << maxi << "," << maxj << ")" << std::endl;
    std::cout << "First Values" << std::endl;

    size_t cnt = 0;
    for(it = m_tank.begin(); it != m_tank.end() && cnt < 5; ++it) {
        std::cout << "(" << it->first.first << "," << it->first.second << ") " << it->second << std::endl;
        cnt++;
    }
}

// =======
// = IOs =
// =======

void sparse_matrice::loadBin( const char *filename )
{
    Maths::ifstream ifs(filename);
    ifs >> Maths::format("old_binary") >> *this;
}

void sparse_matrice::saveBin( const char *filename ) const
{
    Maths::ofstream ofs(filename);
    ofs << Maths::format("old_binary") << *this;
}

void sparse_matrice::loadTxt( const char *filename )
{
    Maths::ifstream ifs(filename);
    ifs >> Maths::format("ascii") >> *this;
}

void sparse_matrice::saveTxt( const char *filename ) const
{
    Maths::ofstream ofs(filename);
    ofs << Maths::format("ascii") << *this;
}

void sparse_matrice::loadMat(const char *filename)
{
    Maths::ifstream ifs(filename);
    ifs >> Maths::format("matlab") >> *this;
}

void sparse_matrice::saveMat( const char *filename ) const
{
    Maths::ofstream ofs(filename);
    ofs << Maths::format("matlab") << *this;
}

void sparse_matrice::load( const char *filename ) {
    try {
        Maths::ifstream ifs(filename);
        ifs >> *this;
    }
    catch (std::string s) {
        std::cout << s << std::endl;
    }
}

void sparse_matrice::save( const char *filename ) const {
    try {
        Maths::ofstream ofs(filename);
        ofs << *this;
    }
    catch (std::string s) {
        std::cout << s << std::endl;
    }
}

std::ostream& operator<<(std::ostream& f,const sparse_matrice &M) {
    sparse_matrice::const_iterator it;
    for(it = M.tank().begin(); it != M.tank().end(); ++it) {
        std::cout << "(" << it->first.first << "," << it->first.second << ") " << it->second << std::endl;
    }
    return f;
}
