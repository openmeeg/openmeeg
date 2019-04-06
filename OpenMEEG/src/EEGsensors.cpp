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

#include <EEGsensors.h>

#include <algorithm>
#include <ciso646>
#include <numeric>     // std::iota
#include <vector>      // std::vector
#include <tuple>

#include <danielsson.h>

namespace OpenMEEG {

    // EEGSensors --------------------------------------
    void EEGSensors::info(int n_lines) const {
        int nb_to_display = (int)std::min((int)m_nb,(int)n_lines);
        std::cout << "EEG electrodes" << std::endl;
        Sensors::info(nb_to_display);
    }

    void EEGSensors::load(const char* filename) {
        // line with
        // 3 elements = sensor position
        // or
        // 4 elements = label + sensor position

        std::ifstream in;
        in.open(filename, std::ios::in);
        std::string s, buf;
        Strings names;
        bool labeled;
        std::set<std::string> labels;
        size_t nlin;
        size_t ncol;
        // determine number of lines, columns and labeled or not
        std::tie(nlin, ncol, labeled) = pre_parse_stream(in);
        if (ncol == 4) {
            labeled = true;
            ncol--;
        }
        else
            labeled = false;
        in.clear();
        in.seekg(0,std::ios::beg); // move the get pointer to the beginning of the file.
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
        // Sensor index
        std::iota(m_pointSensorIdx.begin(), m_pointSensorIdx.end(), 0);
        m_nb = nlin;
    }

    void EEGSensors::save(const char* filename) {
        std::ofstream outfile(filename);
        for ( size_t i = 0; i < getNumberOfPositions(); ++i) {
            // if it has names
            if (hasLabels())
                outfile << m_labels[m_pointSensorIdx[i]] << " ";
            outfile << m_positions.getlin(i) << std::endl;
        }
        return;
    }
}
