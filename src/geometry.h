/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
lastrevision      : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#ifndef H_geometry
#define H_geometry

#include <assert.h>
#include "mesh3.h"

/** \brief  Geometry
    
    Geometry Class
    
**/

class Geometry
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
