// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the GPL open source license.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <string>

#include "Defs.h"
#include "mesh_from_levelset_function.h"

namespace OpenMEEG::OMCGAL {

    Mesh refine_mesh(const Mesh& m_in,const double radius_bound,const double distance_bound,const std::string& sizing_field);
    Mesh decimate_mesh(const Mesh& m_in,const unsigned nb_points);

    Mesh mesh_from_labeled_image(const std::string& input_filename,const double radius_bound,const double distance_bound);
    Mesh mesh_from_levelset_image(const std::string& input_filename,const double levelset_value,const bool positive_inside,const double radius_bound,const double distance_bound);
}
