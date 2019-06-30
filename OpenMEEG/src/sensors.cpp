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

#include <ciso646>
#include <iterator>    // std::distance
#include <numeric>     // std::iota
#include <vector>      // std::vector

#include <danielsson.h>

namespace OpenMEEG {

    bool Sensors::hasSensor(std::string name) const {
        return (std::find(m_labels.cbegin(), m_labels.cend(), name) != m_labels.cend());
    }

    size_t Sensors::getSensorIdx(std::string name) const {
        auto it = std::find(m_labels.cbegin(), m_labels.cend(), name);
        if (it == m_labels.cend()) {
            std::cerr << "Unknown sensor : " << name << std::endl;
            exit(1);
        }
        return std::distance(m_labels.cbegin(), it);
    }

    void Sensors::info(int n_lines) const {
        std::cout << m_type << " sensors:" << std::endl;
        size_t nb_to_display = (int)std::min((int)m_nb,(int)n_lines);
        std::cout << "Nb of sensors : " << m_nb << std::endl;
        if (hasLabels()) {
            std::cout << "Labels" << std::endl;
            for(size_t i = 0; i < nb_to_display; ++i) {
                std::cout << m_labels[i] << std::endl;
            }
            if(m_nb > nb_to_display) {
                std::cout << "..." << std::endl;
            }
        }
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
    }
}
