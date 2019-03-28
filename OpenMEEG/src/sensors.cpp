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

#include <sensors.h>

#include <algorithm>
#include <ciso646>
#include <iterator>     // std::distance
#include <vector>       // std::vector

#include <danielsson.h>

namespace OpenMEEG {

    bool Sensors::hasSensor(std::string name) const {
        return (std::find(m_names.cbegin(), m_names.cend(), name) != m_names.cend());
    }

    size_t Sensors::getSensorIdx(std::string name) const {
        auto it = std::find(m_names.cbegin(), m_names.cend(), name);
        if (it == m_names.cend()) {
            std::cerr << "Unknown sensor : " << name << std::endl;
            exit(1);
        }
        return std::distance(m_names.cbegin(), it);
    }

    void Sensors::load(const char* filename, char filetype) {
        std::ifstream in;
        if(filetype == 't') {
            in.open(filename,std::ios::in);
        } else if(filetype == 'b') {
            in.open(filename,std::ios::in|std::ios::binary);
        } else {
            std::cerr << "ERROR: unkown filetype. " << std::endl; exit(1);
        }

        if ( !in.is_open() ) {
            std::cerr<<"Error Reading File : " << filename << std::endl;
            exit(1);
        }
        Sensors::load(in);
        in.close();
    }

    void Sensors::load(std::istream &in) {

        in >> io_utils::skip_comments('#');

        std::string s, buf;
        Strings names;
        Strings tokens;
        bool labeled = true;
        size_t nlin = 0;
        size_t ncol = 0;
        size_t i = 0;
        // determine number of lines, columns and labeled or not
        while ( std::getline(in,s) ) {
            if ( !s.empty() ) {
                // Tokenize the line.
                std::stringstream iss(s);
                tokens.clear();
                while (iss >> buf) {
                    tokens.push_back(buf);
                    if ( i==0 ) {
                        ncol++;
                    }
                }
                // it is labeled unless there exists a float (i.e containing one '.')
                if (std::count(tokens[0].cbegin(), tokens[0].cend(), '.') == 1) {
                    labeled = false;
                }
                if ( tokens.size() != ncol ) {
                    std::cout << tokens.size() << " != " << ncol << std::endl;
                    std::cerr << "Problem while reading Sensors file" << std::endl;
                    std::cerr << "Each line should have the same number of elements" << std::endl;
                    exit(1);
                }
                ++nlin;
                ++i;
            }
        }
        in.clear();
        in.seekg(0,std::ios::beg); // move the get pointer to the beginning of the file.
        in >> io_utils::skip_comments('#');

        if ( labeled ) {
            ncol--;
        }

        Matrix mat(nlin, ncol);
        i = 0;
        while ( std::getline(in,s) ) {
            if ( !s.empty() ) {
                // Tokenize the line.
                std::stringstream iss(s);
                if ( labeled ) {
                    iss >> buf;
                    names.push_back(buf);
                }
                Vector v(ncol);
                iss >> v;
                mat.setlin(i, v);
                ++i;
            }
        }

        // init private members :
        // positions
        m_positions = mat.submat(0,nlin,0,3);
        // weights
        if (m_geo) { // EIT
            if (ncol == 4) { // if radii were specified
                m_radii = mat.getcol(mat.ncol()-1);
            } else {
                m_radii = Vector(nlin);
                m_radii.set(0.);
            }
            // find triangles on which to inject the currents and compute weights
            findInjectionTriangles();
        } else if ((ncol == 4) && (m_geo == NULL) ) {
            std::cerr << "Sensors:: please specify at constructor stage the geometry on which to apply the spatially extended EIT sensors." << std::endl;
            exit(1);
        } else if (ncol == 7) { // MEG
            m_weights = mat.getcol(mat.ncol()-1);
        } else { // Others
            m_weights = Vector(nlin);
            m_weights.set(1.);
        }
        m_pointSensorIdx = std::vector<size_t>(nlin);

        // orientations
        if ( ncol >= 6 ) {
            m_orientations = mat.submat(0,nlin,3,3);
        }

        // Sensor index
        m_nb = 0;
        if ( labeled ) {
            for ( i = 0; i < nlin; ++i) {
                if ( hasSensor(names[i]) ) {
                    m_pointSensorIdx[i] = getSensorIdx(names[i]);
                } else {
                    m_names.push_back(names[i]);
                    m_pointSensorIdx[i] = m_nb;
                    m_nb++;
                }
            }
        } else {
            for ( i = 0; i < nlin; ++i) {
                m_pointSensorIdx[i] = m_nb;
                m_nb++;
            }
        }
    }

    void Sensors::save(const char* filename) {
        std::ofstream outfile(filename);
        for(size_t i = 0; i < getNumberOfPositions(); ++i) {
            // if it has names
            if (hasNames())
                outfile << m_names[m_pointSensorIdx[i]] << " ";
            outfile << m_positions.getlin(i) << " ";
            // if it has orientations
            if (hasOrientations())
                outfile << m_orientations.getlin(i) << " ";
            // if it has weights (other than 1)
            if (not almost_equal(m_weights.sum(), static_cast<double>(m_weights.size()))) {
                outfile << m_weights(i) << std::endl;
            } else {
                outfile << std::endl;
            }
        }
        return;
    }

    SparseMatrix Sensors::getWeightsMatrix() const {
        SparseMatrix weight_matrix(getNumberOfSensors(),getNumberOfPositions());
        for(size_t i = 0; i < getNumberOfPositions(); ++i) {
            weight_matrix(m_pointSensorIdx[i],i) = m_weights(i);
        }
        return weight_matrix;
    }

    void Sensors::findInjectionTriangles() {
        om_error(m_geo!=NULL);
        m_weights = Vector(m_positions.nlin());
        m_weights.set(1.);
        //To count the number of points that have been mapped to each mesh.
        Strings ci_mesh_names;
        std::vector<size_t>      ci_triangles;

        for ( size_t idx = 0; idx < m_positions.nlin(); ++idx) {
            Triangles triangles;
            const Vect3 current_position(m_positions(idx, 0), m_positions(idx, 1), m_positions(idx, 2));
            Vect3 current_alphas; //not used here
            Triangle current_nearest_triangle; // to hold the closest triangle to electrode.

            double dist;
            std::string s_map=dist_point_geom(current_position, *m_geo, current_alphas, current_nearest_triangle, dist);
            Strings::iterator sit=std::find(ci_mesh_names.begin(),ci_mesh_names.end(),s_map);
            if(sit!=ci_mesh_names.end()){
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

    void Sensors::info() const {
        size_t nb_to_display = (int)std::min((int)m_nb,(int)5);
        std::cout << "Nb of sensors : " << m_nb << std::endl;
        std::cout << "Positions" << std::endl;
        for(size_t i = 0; i < nb_to_display ; ++i) {
            for (size_t j=0;j<m_positions.ncol();++j) {
                std::cout << m_positions(i,j) << " ";
            }
            std::cout << std::endl;
        }
        if(m_nb > nb_to_display) {
            std::cout << "..." << std::endl;
        }

        if(hasOrientations()) {
            std::cout << "Orientations" << std::endl;
            for(size_t i = 0; i < nb_to_display ; ++i) {
                for (size_t j=0;j<m_orientations.ncol();++j) {
                    std::cout << m_orientations(i,j) << " ";
                }
                std::cout << std::endl;
            }
            if(m_nb > nb_to_display) {
                std::cout << "..." << std::endl;
            }
        }
        if(hasRadii()) {
            std::cout << "Radii" << std::endl;
            for(size_t i = 0; i < nb_to_display ; ++i) {
                std::cout << m_radii(i) << " " << std::endl;
            }
            if(m_nb > nb_to_display) {
                std::cout << "..." << std::endl;
            }
        }
        if(hasNames()) {
            std::cout << "Names" << std::endl;
            for(size_t i = 0; i < nb_to_display; ++i) {
                std::cout << m_names[i] << std::endl;
            }
            if(m_nb > nb_to_display) {
                std::cout << "..." << std::endl;
            }
        }
    }
}
