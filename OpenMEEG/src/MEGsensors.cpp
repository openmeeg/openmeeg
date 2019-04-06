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

#include <MEGsensors.h>

#include <algorithm>
#include <ciso646>
#include <numeric>     // std::iota
#include <vector>      // std::vector
#include <tuple>

#include <danielsson.h>

namespace OpenMEEG {
    // MEGSensors --------------------------------------

    void MEGSensors::info(int n_lines) const {
        int nb_to_display = (int)std::min((int)m_nb,(int)n_lines);
        std::cout << "MEG electrodes" << std::endl;
        Sensors::info(nb_to_display);

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
    }

    SparseMatrix MEGSensors::getWeightsMatrix() const {
        SparseMatrix weight_matrix(getNumberOfSensors(),getNumberOfPositions());
        for(size_t i = 0; i < getNumberOfPositions(); ++i) {
            weight_matrix(m_pointSensorIdx[i],i) = m_weights(i);
        }
        return weight_matrix;
    }

    Vector MEGSensors::getOrientation(size_t idx) const {
        return m_orientations.getlin(idx);
    }

    void MEGSensors::setOrientation(size_t idx, Vector& orient) {
        return m_orientations.setlin(idx,orient);
    }

    void MEGSensors::load(const char* filename) {
        // line with
        // 3 elements = position
        // or
        // 4 elements = label + position
        // or
        // 6 elements = position + orientation
        // or
        // 7 elements = label + position + orientation
        // or
        // 8 elements = label + position + orientation + weight

        std::ifstream in;
        in.open(filename, std::ios::in);
        std::string s, buf;
        Strings labels;
        bool labeled;
        size_t nlin;
        size_t ncol;
        // determine number of lines, columns and labeled or not
        std::tie(nlin, ncol, labeled) = pre_parse_stream(in);
        if ((ncol == 4) or (ncol == 7) or (ncol == 8)) {
            labeled = true;
            ncol--;
        }
        else
            labeled = false;
        in >> io_utils::skip_comments('#');

        Matrix mat(nlin, ncol);
        size_t i = 0;
        while ( std::getline(in, s) ) {
            if ( !s.empty() ) {
                // Tokenize the line.
                std::stringstream iss(s);
                if ( labeled ) {
                    iss >> buf;
                    labels.push_back(buf);
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
        // orientations
        if ( ncol >= 6 ) {
            m_orientations = mat.submat(0,nlin,3,3);
        }
        if (ncol == 7) {
            m_weights = mat.getcol(mat.ncol()-1);
        } else {
            m_weights = Vector(nlin);
            m_weights.set(1.);
        }

        m_pointSensorIdx.resize(nlin);
        // Sensor index
        m_nb = 0;
        if ( labeled ) {
            for ( i = 0; i < nlin; ++i) {
                if ( hasSensor(labels[i]) ) {
                    m_pointSensorIdx[i] = getSensorIdx(labels[i]);
                } else {
                    m_labels.push_back(labels[i]);
                    m_pointSensorIdx[i] = m_nb;
                    m_nb++;
                }
            }
        } else {
            std::iota(m_pointSensorIdx.begin(), m_pointSensorIdx.end(), 0);
            m_nb = nlin;
        }
    }

    void MEGSensors::save(const char* filename) {
        std::ofstream outfile(filename);
        for ( size_t i = 0; i < getNumberOfPositions(); ++i) {
            // if it has names
            if (hasLabels())
                outfile << m_labels[m_pointSensorIdx[i]] << " ";
            outfile << m_positions.getlin(i) << " ";
            // if it has orientations
            if (hasOrientations())
                outfile << m_orientations.getlin(i) << " ";
            outfile << std::endl;
        }
        return;
    }
}
