/*! \file
\brief file containing some utility functions.
*/
#ifndef _utils_h
#define _utils_h

#include "setup.h"

#if WIN32
//! allows the use of M_PI (math.h).
#define _USE_MATH_DEFINES
#endif

// #define M_PI       3.14159265358979323846

#include <cstring>
#include <math.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <complex>
//! defines complexe from complex (STL) and reel (setup.h).
typedef std::complex<reel> complexe;

#include <assert.h>

#if WIN32
//! Prevent windows.h from defining min and max
#define NOMINMAX
#include <windows.h>
/*! create the windows version of usleep(int us).
\param us is the sleeping time in microseconds.
*/
inline void usleep(int us) { Sleep(us/1000); }
#else
#include <unistd.h>
#endif

//! start the timer.
void timer_start();
/*! read the timer.
    \attention timer_start() has to be called before.
    \return number of seconds since the call of timer_start()
*/
reel timer_lap();
/*! Initializes the random generator with some given seed
    \attention seed must not be -1
*/
void init_random(int seed);
/*! \return a random value in [0,1]
    \attention Automatically initializes the generator if init_random() has not been called
*/
reel drandom(); // 0 a 1
/*!	\return a random value issued from a gaussian law which parameters are:
    - mean = 0
    - variance = 1
    \attention Automatically initializes the generator if init_random() has not been called
*/
reel gaussienne();
/*! \param name input string containing the name of the file.
    \param extension output string in which the extension is written (memory has to be preallocated).
    \return 0: no error \n
    1: error
*/
int getNameExtension ( const char* name, char* extension);
/*! Find the basename and the extension of a file
    \ingroup loadNd
    \param name image name.
    \param suf extension (output variable, memory has to be allocated).
    \return basename of the file or 0 if no extension.
    \attention The output is static.
*/
char* basename(const char* name,char *& suf);
/*!	\return sign(x) (0, -1 or +1)
*/
inline reel sign(reel x) { return ((x==0)? (reel)0:((x>0)?(reel)1:(reel)-1)); }
#endif

#if WIN32
#include <time.h>
#else
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#endif

static reel ref;
inline void timer_start()
{
#if WIN32
    ref=reel(clock())/CLOCKS_PER_SEC;
#else
    struct tms b;
    times(&b);
    ref=reel(b.tms_utime)/sysconf(_SC_CLK_TCK);
#endif
}

inline reel timer_lap()
{
#if WIN32
    reel lap=reel(clock())/CLOCKS_PER_SEC-ref;
#else
    struct tms b;
    times(&b);
    reel lap=reel(b.tms_utime)/sysconf(_SC_CLK_TCK)-ref;
#endif
    return lap;
}

inline void init_random(int seed) {
    static bool first=true;
    if (seed==-1 && !first)
        return;
    first=false;
    // srand((unsigned int)((seed==-1)?time(0):seed));
    srand(0);
    rand(); // Le premier est biaisé!
}

inline reel drandom()
{
    init_random(-1);
    return reel(rand())/RAND_MAX;
}

inline reel gaussienne()
{
    reel x; // 1er tirage (rejet si nul)
    do 
        x=drandom();
    while (x==0);
    return (reel)(sqrt(-2*log(x))*cos(2*M_PI*drandom()));
}

inline int getNameExtension ( const char* name, char* extension )
{
    const char *point=strrchr(name,'.');
    strcpy(extension,point+1);

    if(point) return 0; else return 1;
};

inline char* basename(const char* name,char *& suf)
{
    static char n[1024];
    strcpy(n,name);
    for (suf=n+strlen(n)-1;*suf!='.' && suf>n;suf--)
        ;
    if (suf==n)
        return 0; // No suffix?
    *suf=0;
    suf++;
    return n;
}

