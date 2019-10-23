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

#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include <IOUtils.H>
#include <vector.h>
#include <matrix.h>
#include <geometry.h>
#include <om_common.h>

#include <OpenMEEG_Export.h>

namespace OpenMEEG {

    /*!
     *  Sensors class for EEG and MEG sensors.
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
        Sensors(): m_nb(0), m_geo(NULL) {} /*!< Default constructor. Number of sensors = 0. */
        Sensors(const Geometry& g): m_nb(0), m_geo(&g) {} /*!< Default constructor with a geometry. Number of sensors = 0. */
        Sensors(const char* filename): m_geo(NULL) { this->load(filename,'t'); } /*!< Construct from file. Option 't' is for text file.*/
        Sensors(const char* filename, const Geometry& g): m_geo(&g) { this->load(filename,'t'); }; /*!< Construct from file and geometry (for EIT). */

        Sensors(const Strings &labels, const Matrix& positions, const Matrix& orientations, const Vector &weights, const Vector &radii);
        Sensors(const Strings &labels, const Matrix& positions, const Matrix& orientations, const Vector &weights, const Vector &radii, const Geometry &g);

        void load(const char* filename, char filetype = 't' ); /*!< Load sensors from file. Filetype is 't' for text file or 'b' for binary file. */
        void load(std::istream &in); /*!< Load description file of sensors from stream. */
        void save(const char* filename);

        size_t getNumberOfSensors() const { return m_nb; } /*!< Return the number of sensors. */
        size_t getNumberOfPositions() const { return m_positions.nlin(); } /*!< Return the number of integration points. */

        Matrix& getPositions() { return m_positions ; } /*!< Return a reference on sensors positions. */
        Matrix getPositions() const { return m_positions ; } /*!< Return a copy of sensors positions */

        Matrix& getOrientations() {return m_orientations ; } /*!< Return a reference on sensors orientations. */
        Matrix getOrientations() const {return m_orientations ; } /*!< Return a copy of sensors orientations. */

        Strings& getNames() {return m_names ; } /*!< Return a reference on sensors names. */
        Strings  getNames() const {return m_names ; } /*!< Return a copy of sensors names. */

        bool hasRadii() const { return m_radii.nlin() > 0 ;} /*!< Return true if contains radii */
        bool hasOrientations() const { return m_orientations.nlin() > 0 ;} /*!< Return true if contains orientations */
        bool hasNames() const { return m_names.size() == m_nb ;} /*!< Return true if contains all sensors names */
        Vector getPosition(size_t idx) const; /*!< Return the position (3D point) of the integration point i. */
        Vector getOrientation(size_t idx) const; /*!< Return the orientations (3D point) of the integration point i. */
        std::string getName(size_t idx) const{ om_assert(idx < m_names.size()); return m_names[idx]; } /*!< Return the name of the idx_th sensor */
        void setPosition(size_t idx, Vector& pos); /*!< Set the position (3D point) of the integration point i. */
        void setOrientation(size_t idx, Vector& orient); /*!< Set the orientation (3D point) of the integration point i. */

        bool hasSensor(std::string name) const;
        size_t getSensorIdx(std::string name) const;
        Triangles getInjectionTriangles(size_t idx) const { om_assert(idx < m_triangles.size()); return m_triangles[idx]; } /*!< For EIT, get triangles under the current injection electrode. */

        Vector getRadii()   const { return m_radii; }
        Vector getWeights() const { return m_weights; }

        SparseMatrix getWeightsMatrix() const;

        bool isEmpty() { if(m_nb == 0) return true; else return false; } /*!< Return if the sensors object is empty. The sensors object is empty if its number of sensors is null. */
        void info() const; /*!< \brief get info about sensors. */

    private:
        size_t m_nb;                        /*!< Number of sensors. */
        Strings m_names;                    /*!< List of sensors names. */
        Matrix m_positions;                 /*!< Matrix of sensors positions. ex: positions(i,j) with  j in {0,1,2} for sensor i */
        Matrix m_orientations;              /*!< Matrix of sensors orientations. ex: orientation(i,j) with  j in {0,1,2} for sensor i */
        Vector m_weights;                   /*!< Weights of integration points */
        Vector m_radii;                     /*!< Areas of the EIT sensors */
        std::vector<Triangles> m_triangles; /*!< Triangles under each EIT sensors */
        const Geometry * m_geo;             /*!< Geometry on which are applied EIT sensors */
        std::vector<size_t> m_pointSensorIdx; /*!< Correspondance between point id and sensor id */
        void findInjectionTriangles();      /*!< Get the triangles under each EIT sensors */
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

    inline Sensors::Sensors(const Strings &labels, const Matrix& positions, const Matrix& orientations, const Vector &weights, const Vector &radii) :
        m_nb(labels.size()), m_names(labels), m_positions(positions), m_orientations(orientations),m_weights(weights), m_radii(radii)
    {
        std::cout << "const" << labels.size() << std::endl;
        m_pointSensorIdx = std::vector<size_t>(labels.size());
        for ( std::size_t i = 0; i < labels.size(); ++i) {
            m_pointSensorIdx[i] = getSensorIdx(m_names[i]);
        }
    }

    inline Sensors::Sensors(const Strings &labels, const Matrix& positions, const Matrix& orientations, const Vector &weights, const Vector &radii, const Geometry &g) :
        m_nb(labels.size()), m_names(labels), m_positions(positions), m_orientations(orientations),m_weights(weights), m_radii(radii), m_geo(&g)
    {
        // find triangles on which to inject the currents and compute weights
        findInjectionTriangles();

        m_pointSensorIdx = std::vector<size_t>(labels.size());
        for ( std::size_t i = 0; i < labels.size(); ++i) {
            m_pointSensorIdx[i] = getSensorIdx(m_names[i]);
        }
    }

}
