// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <string>
#include <filesystem>
#include <algorithm>
#include <cctype>

namespace OpenMEEG {

    inline std::string
    getFilenameExtension(const std::string& name) {
        const std::filesystem::path p(name);
        const std::string ext = p.extension();
        if (ext=="")
            return "";
        return ext.substr(1);
    }

    inline std::string
    tolower(const std::string& s) {
        std::string res = s;
        std::transform(res.begin(),res.end(),res.begin(),
                       [](unsigned char c){ return static_cast<unsigned char>(std::tolower(c)); });
        return res;
    }
}
