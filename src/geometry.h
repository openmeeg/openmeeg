#ifndef H_geometry
#define H_geometry

#include <assert.h>
#include "mesh3.h"

class geometry
{
    int n;
    double *sigin;
    double *sigout;
    mesh *M;

public:
    geometry() {n=0;}
    ~geometry() {destroy();}
    int read(char* geomFileName, char* condFileName);
    inline int nb() const {return n;}
    inline double sigma_in(int i) const {return (i<n)?sigin[i]:0;} // 0 for sig(n)
    inline double sigma_out(int i) const {return (i<n)?sigout[i]:0;} // 0 for sig(n)
    inline double sigma(int i) const {return (i<n)?sigin[i]:0;} // 0 for sig(n)
    inline const mesh &surf(int i) const {return M[i];}
    inline mesh& getM(int i) {assert(i>-1 && i<n); return M[i];}

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
