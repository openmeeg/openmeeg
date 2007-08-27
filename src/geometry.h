#ifndef H_geometry
#define H_geometry

#include <assert.h>
#include "mesh3.h"

/** \brief  Geometry
    
    Geometry Class
    
© Copyright 2007-2007 Odyssee INRIA . All Rights Reserved.

    \author $LastChangedBy$
    \date $LastChangedDate$
    \version $Rev$  \sa
**/

class Geometry
{
    int n;
    double *sigin;
    double *sigout;
    Mesh *M;

public:
    Geometry() {n=0;}
    ~Geometry() {destroy();}
    int read(char* geomFileName, char* condFileName);
    inline int nb() const {return n;}
    inline double sigma_in(int i) const {return (i<n)?sigin[i]:0;} // 0 for sig(n)
    inline double sigma_out(int i) const {return (i<n)?sigout[i]:0;} // 0 for sig(n)
    inline double sigma(int i) const {return (i<n)?sigin[i]:0;} // 0 for sig(n)
    inline const Mesh &surf(int i) const {return M[i];}

    inline       Mesh& getM(const int i)       {assert(i>-1 && i<n); return M[i];}
    inline const Mesh& getM(const int i) const {assert(i>-1 && i<n); return M[i];}

private:
    void destroy() {
        if (n!=0){
            n=0;
            delete []M;
            delete []sigin;
            delete []sigout;
        }
    }

};

#endif
