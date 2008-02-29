/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#ifndef MESHDESCRIPTION_NORMALS_H
#define MESHDESCRIPTION_NORMALS_H

namespace MeshDescription {

    struct NoNormals {
        typedef NoNormals PixelType;

        template <typename T>
        const PixelType& operator()(const T&) const { return *this; }
    };

    inline std::istream&
    operator>>(std::istream& is,NoNormals&) {
        #ifdef DEBUG
        std::cerr << "Normals are not compiled in: not loading normals." << std::endl;
        #endif
        return is;
    }
}

#endif  //  ! MESHDESCRIPTION_NORMALS_H
