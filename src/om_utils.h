#if WIN32
//! allows the use of M_PI (math.h).
#define _USE_MATH_DEFINES
#endif
// #define M_PI       3.14159265358979323846
#include <math.h>

inline int getNameExtension ( const char* name, char* extension )
{
    const char *point=strrchr(name,'.');
    strcpy(extension,point+1);

    if(point) return 0; else return 1;
};

inline void init_random(int seed) {
    static bool first=true;
    if (seed==-1 && !first)
        return;
    first=false;
    // srand((unsigned int)((seed==-1)?time(0):seed));
    srand(0);
    rand(); // the first is biased!
}

inline double drandom()
{
    init_random(-1);
    return double(rand())/RAND_MAX;
}

inline double gaussienne()
{
    double x;
    do 
        x=drandom();
    while (x==0);
    return (double)(sqrt(-2*log(x))*cos(2*M_PI*drandom()));
}
