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

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <vector>

#include <IOUtils.H>
#include <vector.h>
#include <matrix.h>
#include <om_common.h>

#include <OpenMEEG_Export.h>

namespace OpenMEEG {

    /*!
     *  Sensors abstract class for EEG MEG, EIT, ECoG sensors.
     *  This class is made for reading sensors description file. This description file is a file text. Sensors may have names (labels)
     *  in the first column of the file (it has to contains at least one character to be considered as label)
     *  the file can have the shape of (neglecting if present the first, label column):
     *  <ul>
     *
     *  <li> 1 line per sensor and 3 columns (EEG sensors OR MEG sensors without orientation OR EIT punctual patches)
     *        <ul TYPE="circle">
     *        <li> the 1st, 2nd and 3rd columns are respectively position coordinates x, y, z of sensor  </li>
     *        </ul>
     *  </li>
     *  <li> 1 line per sensor and 4 columns (EEG EIT patches (circular patches)) :
     *        <ul TYPE="circle">
     *        <li> the 1st, 2nd and 3rd are respectively position coordinates x, y, z of sensor  </li>
     *        <li> the 4th is the patche radius (unit relative to the mesh)  </li>
     *        </ul>
     *  </li>
     *  <li> 1 line per sensor and 6 columns (MEG sensors) :
     *        <ul TYPE="circle">
     *        <li> the 1st, 2nd and 3rd are respectively position coordinates x, y, z of sensor  </li>
     *        <li> the 4th, 5th and 6th are coordinates of vector orientation </li>
     *        </ul>
     *  </li>
     *  <li> 1 line per integration point for each sensor and 7 columns (MEG sensors) :
     *        <ul TYPE="circle">
     *        <li> the 1st, 2nd and 3rd are respectively position coordinates x, y, z of sensor  </li>
     *        <li> the 4th, 5th and 6th are coordinates of vector orientation </li>
     *        <li> the 7th is the weight to apply for numerical integration (uses sensor name) </li>
     *        </ul>
     *  </li>
     *  </ul>
     */

    class OPENMEEG_EXPORT Sensors {

    public:
        Sensors(): m_nb(0) {} /*!< Default constructor. Number of sensors = 0. */

        virtual void load(const char* filename) = 0; /*!< Load sensors from file. Filetype is 't' for text file or 'b' for binary file. */
        virtual void save(const char* filename) = 0;

        size_t getNumberOfSensors() const { return m_nb; } /*!< Return the number of sensors. */
        size_t getNumberOfPositions() const { return m_positions.nlin(); } /*!< Return the number of integration points. */

        Matrix& getPositions() { return m_positions ; } /*!< Return a reference on sensors positions. */
        Matrix getPositions() const { return m_positions ; } /*!< Return a copy of sensors positions */

        Strings& getNames() {return m_labels ; } /*!< Return a reference on sensors names. */
        Strings  getNames() const {return m_labels ; } /*!< Return a copy of sensors names. */

        bool hasLabels() const { return m_labels.size() == m_nb ;} /*!< Return true if contains all sensors names */
        Vector getPosition(size_t idx) const; /*!< Return the position (3D point) of the integration point i. */
        std::string getName(size_t idx) const{ om_assert(idx < m_labels.size()); return m_labels[idx]; } /*!< Return the name of the idx_th sensor */
        void setPosition(size_t idx, Vector& pos); /*!< Set the position (3D point) of the integration point i. */

        bool hasSensor(std::string name) const;
        size_t getSensorIdx(std::string name) const;
        bool isEmpty() const { return (m_nb == 0); } /*!< Return if the sensors object is empty. The sensors object is empty if its number of sensors is null. */
        void info(int n_lines = 5) const; /*!< \brief get n_lines first lines info about sensors. */

    protected:
        size_t m_nb;                          /*!< Number of sensors. */
        Strings m_labels;                     /*!< List of sensors names. */
        Matrix m_positions;                   /*!< Matrix of sensors positions. ex: positions(i,j) with  j in {0,1,2} for sensor i */
        std::vector<size_t> m_pointSensorIdx; /*!< Correspondance between point id and sensor id */
    };

    inline Vector Sensors::getPosition(size_t idx) const {
        return m_positions.getlin(idx);
    }

    inline void Sensors::setPosition(size_t idx, Vector& pos) {
        return m_positions.setlin(idx,pos);
    }
}

namespace {
    std::tuple<size_t, size_t, bool> pre_parse_stream(std::istream& in) {
        in >> io_utils::skip_comments('#');

        OpenMEEG::Strings tokens;
        bool seems_labeled = true;
        std::string s, buf;
        size_t nlin = 0, ncol = 0;

        while ( std::getline(in,s) ) {
            if ( !s.empty() ) {
                // Tokenize the line.
                std::stringstream iss(s);
                tokens.clear();
                while (iss >> buf) {
                    tokens.push_back(buf);
                }
                // it is labeled unless there exists a float (i.e containing one '.')
                if (std::count(tokens[0].cbegin(), tokens[0].cend(), '.') == 1) {
                    seems_labeled = false;
                }
                if (nlin == 0) // first line
                    ncol = tokens.size();
                else if ( tokens.size() != ncol ) {
                    std::cerr << tokens.size() << " != " << ncol << std::endl;
                    std::cerr << "Problem while reading sensors file" << std::endl;
                    std::cerr << "Each line should have the same number of elements" << std::endl;
                    exit(1);
                }
                ++nlin;
            }
        }
        if (nlin == 0) {
            std::cerr << "Problem while reading sensors file" << std::endl;
            std::cerr << "File does not contain any sensor" << std::endl;
            exit(1);
        }
        in.clear();
        in.seekg(0,std::ios::beg); // move the get pointer to the beginning of the file.
        return std::make_tuple(nlin, ncol, seems_labeled);
    }
}
