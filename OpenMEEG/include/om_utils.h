/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

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

#pragma once

#if WIN32
#define _USE_MATH_DEFINES
#endif

#include <string>
#include <cmath>
#include <random>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#ifdef USE_OMP
#include <omp.h>
#endif

#include "OpenMEEGConfigure.h"

#ifdef USE_PROGRESSBAR
    #define PROGRESSBAR(a,b) progressbar((a),(b))
#else
    #define PROGRESSBAR(a,b)
#endif

namespace OpenMEEG {

#ifndef M_PI
    constexpr double Pi = 3.14159265358979323846;
#else
    constexpr double Pi = M_PI;
#endif

#if WIN32
    constexpr char PathSeparator[] = "/\\";
#else
    constexpr char PathSeparator = '/';
#endif

    constexpr double MagFactor = 1e-7;

    inline std::string
    getFilenameExtension(const std::string& name) {
        const std::string::size_type idx = name.find_last_of('.');
        if (idx==std::string::npos)
            return "";
        return name.substr(idx+1);
    }

    inline std::string
    tolower(const std::string& s) {
        std::string res = s;
        std::transform(res.begin(),res.end(),res.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        return res;
    }

    /// \return absolute path of file \param name.
    
    inline std::string
    absolute_path(const std::string& name) {
        const std::string::size_type pos = name.find_last_of(PathSeparator);
        return (pos==std::string::npos) ? "" : name.substr(0,pos+1);
    }

    /// \return true if name is a relative path. \param name

    inline bool
    is_relative_path(const std::string& name) {
    #if WIN32
        const char c0 = name[0];
        if (c0=='/' || c0=='\\')
            return false;

        const char c1 = name[1];
        const char c2 = name[2];
        return !(std::isalpha(c0) && c1==':' && (c2=='/' || c2=='\\'));
    #else
        return name[0]!=PathSeparator;
    #endif
    }

    inline void
    disp_argv(const int argc,char **argv) {
        std::cout << std::endl << "| ------ " << argv[0] << std::endl;
        for (int i=1;i<argc;++i)
            std::cout << "| " << argv[i] << std::endl;
        std::cout << "| -----------------------" << std::endl;
    }

#ifdef USE_PROGRESSBAR
    inline void progressbar(const unsigned n,const unsigned N,const unsigned w=20) {
        // w : nb of steps
        const char cprog  = '.';
        const char cprog1 = '*';
        const char cbeg   = '[';
        const char cend   = ']';
        const unsigned p = std::min(static_cast<unsigned>(floor(1.f*n*(w+1)/N)),w);

        static unsigned pprev = -1;
        if (N>1) {
            if (n==0)
                pprev = -1;

            if (p!=pprev) {
                if (n>1) {
                    // clear previous string
                    for (unsigned i=0;i<(w+2);++i)
                        std::cout << "\b";

                    std::cout << cbeg;
                    for (unsigned i=0;i<p;++i)
                        std::cout << cprog1;

                    for (unsigned i=p;i<w;++i)
                        std::cout << cprog;

                    std::cout<< cend;
                }
            }
            pprev = p;
            if (n>=(N-1))
                std::cout << std::endl;
            std::cout.flush();
        }
    }
#endif

    inline void warning(std::string message) {
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        std::cout << "!!!!!!!!!!! WARNING !!!!!!!!!!!" << std::endl;
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        std::cout << message << std::endl;
    }

    inline void print_version(const char* cmd) {
        #ifdef USE_OMP
            std::string omp_support = " using OpenMP\n Executing using " + std::to_string(omp_get_max_threads()) + " threads.";
        #else
            std::string omp_support = "";
        #endif

        std::ostringstream display_info;
        display_info << cmd;
        display_info << " version " << version;
        display_info << " compiled at " << __DATE__ << " " << __TIME__;
        display_info << omp_support;
        std::cout << display_info.str() << std::endl << std::endl;
    }
}
