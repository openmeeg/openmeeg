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
#include <chrono>
#include <sstream>
#include <algorithm>
#include <cctype>

namespace OpenMEEG {

#if WIN32
    constexpr char PathSeparator[] = "/\\";
#else
    constexpr char PathSeparator = '/';
#endif

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

    /// \return absolute path of file \param name .
    
    inline std::string
    absolute_path(const std::string& name) {
        const std::string::size_type pos = name.find_last_of(PathSeparator);
        return (pos==std::string::npos) ? "" : name.substr(0,pos+1);
    }

    /// \return the base name of file \param name .
    
    inline std::string
    basename(const std::string& name) {
        const std::string::size_type pos = name.find_last_of(PathSeparator);
        return (pos==std::string::npos) ? name : name.substr(pos+1);
    }

    /// \return true if name is a relative path. \param name

    inline bool
    is_relative_path(const std::string& name) {
    #if WIN32
        const std::string& sep = PathSeparator;
        const char c0 = name[0];
        if (sep.find(c0)!=std::string::npos)
            return false;
        return !(std::isalpha(c0) && name[1]==':' && sep.find(name[2])!=std::string::npos);
    #else
        return name[0]!=PathSeparator;
    #endif
    }
}
