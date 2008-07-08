/* FILE: $Id: sensors.h 222 2008-04-08 06:14:41Z gramfort $ */

/*
Project Name : OpenMEEG

version           : $Revision: 222 $
last revision     : $Date: 2008-04-08 08:14:41 +0200 (Tue, 08 Apr 2008) $
modified by       : $LastChangedBy: gramfort $
last modified     : $LastChangedDate: 2008-04-08 08:14:41 +0200 (Tue, 08 Apr 2008) $

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

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

#ifndef SENSORS_H
#define SENSORS_H

#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "IOUtils.H"
#include "vector.h"
#include "matrix.h"
#include "symmatrix.h"
#include "sparse_matrix.h"

/*!
 *  Sensors class for EEG and MEG sensors.
 *  This class is made for read sensors description file. This description file is a file text which can take the shape of :
 *  <ul>
 *    <li> 1 line per 1 sensor and 7 columns (MEG sensors) :
 *        <ul TYPE="circle">
 *        <li> the 1st column is sensors names </li>
 *        <li> the 2nd, 3rd and 4th are respectively positions coordinates x, y, z of sensor  </li>
 *        <li> the 5th, 6th and 7th are coordinates of vector orientation </li>
 *        </ul>
 *  </li>
 *  <li> 1 line per 1 sensor and 6 columns (MEG sensors) :
 *        <ul TYPE="circle">
 *        <li>- the 1st, 2nd and 3rd are respectively positions coordinates x, y, z of sensor  </li>
 *        <li>- the 4th, 5th and 6th are coordinates of vector orientation </li>
 *        </ul>
 *  </li>
 *    <li> 1 line per 1 sensor and 4 columns (EEG sensors or MEG sensors without orientation) :
 *        <ul TYPE="circle">
 *        <li>- the 1st column is sensors names </li>
 *        <li>- the 2nd, 3rd and 4th are respectively positions coordinates x, y, z of sensor  </li>
 *        </ul>
 *  </li>
 *    <li> 1 line per 1 sensor and 3 columns (EEG sensors or MEG sensors without orientation) :
 *        <ul TYPE="circle">
 *        <li>- the 1st, 2nd and 3rd are respectively positions coordinates x, y, z of sensor  </li>
 *        </ul>
 *  </li>
 *  </ul>
 */

class Sensors {
private:
    size_t m_nb;                       /*!< Number of sensors. */
    std::vector<std::string> m_names;  /*!< List of sensors names. */
    Matrix m_positions;               /*!< Matrix of sensors positions. ex: positions(i,j) with  j in {0,1,2} for sensor i */
    Matrix m_orientations;            /*!< Matrix of sensors orientations. ex: orientation(i,j) with  j in {0,1,2} for sensor i */
    Vector m_weights;                 /*!< Weights of integration points */
    std::vector<size_t> m_pointSensorIdx; /*!< Correspondance between point id and sensor id */
    void copy(const Sensors& S);       /*!< Copy function. Copy sensor S in current sensor object. ex. senors S1; ...; sensors S2(S1); */

public:
    Sensors(): m_nb(0) {} /*!< Default constructor. Number of sensors = 0. */
    Sensors(char* filename); /*!< Construct from file. Option 't' is for text file, and 'b' is for binary file. */
    Sensors(const Sensors& S) { copy(S); }        /*!< Copy constructor. */
    ~Sensors() { m_nb=0; }                        /*!< Destructor. Number of sensors = 0. */

    Sensors& operator=(const Sensors& S); /*!< Copy operator. Copy sensor S in current sensor object. ex. sensors S1; ...; sensors S2 = S1; */

    void load(const char* filename, char filetype = 't' ); /*!< Load sensors from file. Filetype is 't' for text file or 'b' for binary file. */
    void load(std::istream &in); /*!< Load description file of sensors from stream. */
    void save(const char* filename);

    size_t getNumberOfSensors() const { return m_nb; } /*!< Return the number of sensors. */
    size_t getNumberOfPositions() const { return m_positions.nlin(); } /*!< Return the number of integration points. */

    Matrix& getPositions() { return m_positions ; } /*!< Return a reference on sensors positions. */
    Matrix getPositions() const { return m_positions ; } /*!< Return a copy of sensors positions */

    Matrix& getOrientations() {return m_orientations ; } /*!< Return a reference on sensors orientations. */
    Matrix getOrientations() const {return m_orientations ; } /*!< Return a copy of sensors orientations. */

    std::vector<std::string>& getSensorsNames() {return m_names ; } /*!< Return a reference on sensors ids. */
    std::vector<std::string> getSensorsNames() const {return m_names ; } /*!< Return a copy of sensors ids. */

    bool hasOrientations() const { return m_orientations.nlin() > 0 ;} /*!< Return true if contains orientations */
    bool hasNames() const { return m_names.size() == m_nb ;} /*!< Return true if contains all sensors names */
    Vector getPosition(size_t idx) const; /*!< Return the position (3D point) of the integration point i. */
    Vector getOrientation(size_t idx) const; /*!< Return the orientations (3D point) of the integration point i. */
    void setPosition(size_t idx, Vector& pos); /*!< Set the position (3D point) of the integration point i. */
    void setOrientation(size_t idx, Vector& orient); /*!< Set the orientation (3D point) of the integration point i. */

    bool hasSensor(std::string name);
    size_t getSensorIdx(std::string name);

    SparseMatrix getWeightsMatrix() const;

    bool isEmpty() { if(m_nb == 0) return true; else return false; } /*!< Return if the sensors object is empty. The sensors object is empty if its number of sensors is null. */
};

inline Vector Sensors::getPosition(size_t idx) const {
    return m_positions.getlin(idx);
}

inline Vector Sensors::getOrientation(size_t idx) const {
    return m_orientations.getlin(idx);
}

inline void Sensors::setPosition(size_t idx, Vector& pos) {
    return m_positions.setlin(idx,pos);
}

inline void Sensors::setOrientation(size_t idx, Vector& orient) {
    return m_orientations.setlin(idx,orient);
}

inline std::ostream& operator<<(std::ostream& f,const Sensors &S) {
    f << "Nb of sensors : " << S.getNumberOfSensors() << std::endl;
    f << "Positions" << std::endl;
    f << S.getPositions();
    if(S.hasOrientations()) {
        f << "Orientations" << std::endl;
        f << S.getOrientations();
    }
    if(S.hasNames()) {
        f << "Names" << std::endl;
        std::vector<std::string> names = S.getSensorsNames();
        for(size_t i = 0; i < names.size(); ++i) {
            f << names[i] << std::endl;
        }
    }
    return f;
}


#endif
