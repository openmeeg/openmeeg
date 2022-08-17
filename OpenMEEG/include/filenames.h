// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#ifdef WIN32
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

#ifdef WIN32
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
    #ifdef _MSC_VER
    #pragma warning( push )
    #pragma warning( disable:4244 )  /*  warning C4244: '=': conversion from 'int' to 'char', possible loss of data */
    #endif
    tolower(const std::string& s) {
        std::string res = s;
        std::transform(res.begin(),res.end(),res.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        return res;
    }
    #ifdef _MSC_VER
    #pragma warning( pop )
    #endif

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
        bool is_rel = false;
    #ifdef WIN32
        const std::string& sep = PathSeparator;
        const char c0 = name[0];
        if (sep.find(c0)!=std::string::npos)
            is_rel = false;
        else {
            is_rel = !(std::isalpha(c0) && name[1]==':' && sep.find(name[2])!=std::string::npos);
        }
    #else
        is_rel = name[0]!=PathSeparator;
    #endif
        return is_rel;
    }
}
