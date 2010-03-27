/* FILE: $Id$ */

/*
Project Name : OpenMEEG

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

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <assert.h>
#include "mesh3.h"

/** \brief  Geometry
    
    Geometry Class
    
**/

class OPENMEEG_EXPORT Geometry
{
    int n;
    double *sigin;
    double *sigout;
    Mesh *M;
    size_t m_size; // Number of triangles + Number of points
    bool has_cond;

public:
    Geometry() {n=0;}
    ~Geometry() {destroy();}
    int read(const char* geomFileName, const char* condFileName = NULL);
    inline int nb() const {return n;}
    inline size_t size() const {return m_size;}
    inline double sigma_in(int i) const {return (i<n)?sigin[i]:0;} // 0 for sig(n)
    inline double sigma_out(int i) const {return (i<n)?sigout[i]:0;} // 0 for sig(n)
    inline double sigma(int i) const {return (i<n)?sigin[i]:0;} // 0 for sig(n)
    inline const Mesh &surf(int i) const {return M[i];}

    inline       Mesh& getM(const int i)       {assert(i>-1 && i<n); return M[i];}
    inline const Mesh& getM(const int i) const {assert(i>-1 && i<n); return M[i];}

    inline int getNumberOfPoints() const {
        int number_points_total = 0;
        for(int i = 0; i < n; ++i)
        {
            number_points_total += M[i].nbPts();
        }
        return number_points_total;
    }

    bool selfCheck() const;
    bool check(const Mesh& m) const;
private:
    void destroy() {
        if (n!=0){
            n=0;
            delete []M;
            if(has_cond)
            {
                delete []sigin;
                delete []sigout;
            }
        }
    }

};

#endif
