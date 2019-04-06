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

#include <vector>

#include <IOUtils.H>
#include <vector.h>
#include <matrix.h>
#include <geometry.h>
#include <sensors.h>

#include <OpenMEEG_Export.h>

namespace OpenMEEG {

    /*!
     *  Sensors class for EEG and MEG sensors.
     *  This class is made for reading sensors description file. This description file is a file text. Sensors may have names (labels)
     *  in the first column of the file (it has to contains at least one character to be considered as label)
     *  the file can have the shape of (neglecting if present the first, label column):
     *  <ul>
     *
     *  <li> 1 line per sensor and 3 columns (MEG sensors without orientation)
     *        <ul TYPE="circle">
     *        <li> the 1st, 2nd and 3rd columns are respectively position coordinates x, y, z of sensor  </li>
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

    class OPENMEEG_EXPORT MEGSensors : public Sensors {
    public:
        MEGSensors(): Sensors() { }; /*!< Default constructor. */
        MEGSensors(const char* filename): Sensors() { load(filename); }; /*!< Construct from file. */

        void info(int n_lines = 5) const; /*!< \brief get n_lines first lines info about sensors. */
        void load(const char* filename); /*!< Load sensors from file. */
        void save(const char* filename);
        Matrix& getOrientations() {return m_orientations ; } /*!< Return a reference on sensors orientations. */
        Matrix getOrientations() const {return m_orientations ; } /*!< Return a copy of sensors orientations. */
        bool hasOrientations() const { return m_orientations.nlin() > 0 ;} /*!< Return true if contains orientations */
        Vector getOrientation(size_t idx) const; /*!< Return the orientations (3D point) of the integration point i. */
        void setOrientation(size_t idx, Vector& orient); /*!< Set the orientation (3D point) of the integration point i. */

        Vector getWeights() const { return m_weights; }
        SparseMatrix getWeightsMatrix() const;
    private:
        Vector m_weights;                   /*!< Weights of integration points */
        Matrix m_orientations;              /*!< Matrix of sensors orientations. ex: orientation(i,j) with  j in {0,1,2} for sensor i */
    };
}
