#ifndef _OM_UTILS_H_
#define _OM_UTILS_H_

#if WIN32
#define _USE_MATH_DEFINES
#endif
#include <cstring>
#include <cmath>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

inline void getNameExtension ( const char* name, char* extension )
{
    const char *point = strrchr(name,'.');
    if(point)
    {
        strcpy(extension,point+1);
    } else {
        std::cerr << "Unable to get file extension from file : " << *name << std::endl;
        exit(1);
    }
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

inline double gaussian()
{
    double x;
    do
        x=drandom();
    while (x==0);
    return (double)(sqrt(-2*log(x))*cos(2*M_PI*drandom()));
}

inline void disp_argv(int argc, char **argv) {
    std::cout << std::endl << "| ------ " << argv[0] << std::endl;
    for( int i = 1; i < argc; i += 1 )
    {
        std::cout << "| " << argv[i] << std::endl;
    }
    std::cout << "| -----------------------" << std::endl;
}

#endif /* _OM_UTILS_H_ */

