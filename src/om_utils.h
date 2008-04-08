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

#ifndef _OM_UTILS_H_
#define _OM_UTILS_H_

#if WIN32
#define _USE_MATH_DEFINES
#endif

#include <string>
#include <cmath>
#include <iostream>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define MU0 1.0 //1.25e-6

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

inline void progressbar(int n, int N, int w = 20) {
#ifdef USE_PROGRESSBAR
    // w : nb of steps
    const char* cprog = ".";
    const char* cprog1 = "*";
    const char* cbeg = "[";
    const char* cend = "]";
    int p = std::min( (int)floor(n*(w+1)/N), w);

    static int pprev = -1;
    if (n == 0) {
        pprev = -1;
    }

    if (p != pprev) {
        if (n>1) {
            // clear previous string
            for(size_t i = 0; i < (w+2); ++i)
                printf( "\b" );

            printf( cbeg );
            for(size_t i = 0; i < p; ++i) {
                printf( cprog1 );
            }
            for(size_t i = p; i < w; ++i) {
                printf( cprog );
            }
            printf( cend );
        }
    }
    pprev = p;
    if (n >= (N-1)) {
        printf("\n");
    }
    std::cout.flush();
#endif
}

inline void warning(std::string message) {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "!!!!!!!!!!! WARNING !!!!!!!!!!!" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << message << std::endl;
}

#endif /* _OM_UTILS_H_ */

