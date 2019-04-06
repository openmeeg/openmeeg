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

#include <EITsensors.h>

#include <algorithm>
#include <ciso646>
#include <iterator>    // std::distance
#include <numeric>     // std::iota
#include <vector>      // std::vector
#include <tuple>

#include <danielsson.h>

namespace OpenMEEG {
    // EITSensors --------------------------------------
    void EITSensors::info(int n_lines) const {
        int nb_to_display = (int)std::min((int)m_nb,(int)n_lines);
        std::cout << "EIT electrodes" << std::endl;
        Sensors::info(nb_to_display);

        if (hasRadii()) {
            std::cout << "Radii" << std::endl;
            for(size_t i = 0; i < nb_to_display ; ++i) {
                std::cout << m_radii(i) << " " << std::endl;
            }
            if(m_nb > nb_to_display) {
                std::cout << "..." << std::endl;
            }
        }
    }

    void EITSensors::load(const char* filename) {
        // line with
        // 3 elements = position
        // or
        // 4 elements = label + position
        // or
        // 4 elements = position + radius
        // or
        // 5 elements = label + position + radius

        std::ifstream in;
        in.open(filename, std::ios::in);
        om_error(m_geo!=NULL);
        std::string s, buf;
        Strings names;
        bool labeled;
        std::set<std::string> labels;
        size_t nlin;
        size_t ncol;
        // determine number of lines, columns and labeled or not
        std::tie(nlin, ncol, labeled) = pre_parse_stream(in);
        if ((ncol == 5) or labeled) {
            labeled = true;
            ncol--;
        }
        in >> io_utils::skip_comments('#');

        Matrix mat(nlin, ncol);
        size_t i = 0;
        while ( std::getline(in, s) ) {
            if ( !s.empty() ) {
                // Tokenize the line.
                std::stringstream iss(s);
                if ( labeled ) {
                    iss >> buf;
                    m_labels.push_back(buf);
                    labels.insert(buf);
                }
                Vector v(ncol);
                iss >> v;
                mat.setlin(i, v);
                ++i;
            }
        }
        if ( labeled and labels.size() != nlin) {
            std::cerr << "Problem while reading the EEG sensors file" << std::endl;
            std::cerr << "Each label should be unique" << std::endl;
            exit(1);
        }

        // init private members :
        // positions
        m_positions = mat.submat(0,nlin,0,3);
        m_pointSensorIdx.resize(nlin);
        // radii
        if (not m_geo) {
            std::cerr << "Sensors:: should not happen: geometry is undefined." << std::endl;
            exit(1);
        }
        if (ncol == 4) { // if radii were specified
            m_radii = mat.getcol(mat.ncol()-1);
        } else {
            m_radii = Vector(nlin);
            m_radii.set(0.);
        }
        // find triangles on which to inject the currents and compute weights
        findInjectionTriangles();
        // Sensor index
        std::iota(m_pointSensorIdx.begin(), m_pointSensorIdx.end(), 0);
        m_nb = nlin;
    }

    void EITSensors::save(const char* filename) {
        std::ofstream outfile(filename);
        bool has_radii = not almost_equal(m_radii.sum(), 0.);
        for ( size_t i = 0; i < getNumberOfPositions(); ++i) {
            // if it has names
            if (hasLabels())
                outfile << m_labels[m_pointSensorIdx[i]] << " ";
            outfile << m_positions.getlin(i) << " ";
            // if it has radii (other than 0)
            if (has_radii) {
                outfile << m_radii(i) << std::endl;
            } else {
                outfile << std::endl;
            }
        }
        return;
    }

    void EITSensors::findInjectionTriangles() {
        om_error(m_geo!=NULL);
        m_weights = Vector(m_positions.nlin());
        m_weights.set(1.);
        // To count the number of points that have been mapped to each mesh.
        Strings ci_mesh_names;
        std::vector<size_t> ci_triangles;

        for ( size_t idx = 0; idx < m_positions.nlin(); ++idx) {
            Triangles triangles;
            const Vect3 current_position(m_positions(idx, 0), m_positions(idx, 1), m_positions(idx, 2));
            Vect3 current_alphas; //not used here
            Triangle current_nearest_triangle; // to hold the closest triangle to electrode.

            double dist;
            std::string s_map=dist_point_geom(current_position, *m_geo, current_alphas, current_nearest_triangle, dist);
            Strings::iterator sit=std::find(ci_mesh_names.begin(),ci_mesh_names.end(),s_map);
            if (sit!=ci_mesh_names.end()){
                size_t idx2=std::distance(ci_mesh_names.begin(),sit);
                ci_triangles[idx2]++;
            }
            else{
                ci_mesh_names.push_back(s_map);
                ci_triangles.push_back(1);
            }

            triangles.push_back(current_nearest_triangle);
            std::set<size_t> index_seen; // to avoid infinite looping
            index_seen.insert(current_nearest_triangle.index());
            if ( not almost_equal(m_radii(idx), 0.) ) {
                // if the electrode is larger than the triangle, look for adjacent triangles
                if ( current_nearest_triangle.area() < 4.*M_PI*std::pow(m_radii(idx),2) ) {
                    std::stack<Triangle *> tri_stack;
                    tri_stack.push(&current_nearest_triangle);
                    while ( !tri_stack.empty() ) {
                        Triangle * t = tri_stack.top();
                        tri_stack.pop();
                        if ( (t->center()-current_position).norm() < m_radii(idx) ) {
                            if (t->index() != current_nearest_triangle.index()) //don't push the nearest triangle twice
                                triangles.push_back(*t);
                            Interface::VectPTriangle t_adj = m_geo->interface(s_map).adjacent_triangles(*t);
                            if ( index_seen.insert(t_adj[0]->index()).second ) tri_stack.push(t_adj[0]);
                            if ( index_seen.insert(t_adj[1]->index()).second ) tri_stack.push(t_adj[1]);
                            if ( index_seen.insert(t_adj[2]->index()).second ) tri_stack.push(t_adj[2]);
                        }
                    }
                }
                // now set the weight as the ratio between the wanted sensor surface and the actual surface
                // (should be close to 1)
                double triangles_area = 0.;
                for ( Triangles::const_iterator tit = triangles.begin(); tit != triangles.end(); ++tit) {
                    triangles_area += tit->area();
                }
                m_weights(idx) = M_PI * std::pow(m_radii(idx),2) / triangles_area;
            }
            m_triangles.push_back(triangles);
        }
        for(size_t i=0;i<ci_mesh_names.size();++i)
            std::cout<<ci_triangles[i]<<" points have been mapped to mesh "<<ci_mesh_names[i]<<std::endl;
    }

    SparseMatrix EITSensors::getWeightsMatrix() const {
        SparseMatrix weight_matrix(getNumberOfSensors(),getNumberOfPositions());
        for(size_t i = 0; i < getNumberOfPositions(); ++i) {
            weight_matrix(m_pointSensorIdx[i],i) = m_weights(i);
        }
        return weight_matrix;
    }

}
