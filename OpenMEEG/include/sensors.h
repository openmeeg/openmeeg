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
#include <sparse_matrix.h>

#include <geometry.h>
#include <om_common.h>

#include <OpenMEEG_Export.h>

namespace OpenMEEG {

    /// Sensors class for EEG and MEG sensors.
    /// This class is made for reading sensors description file. This description file is a file text.
    /// Sensors may have names (labels) in the first column of the file (it has to contains at least
    /// one character to be considered as label) the file can have the shape of (neglecting if present
    /// the first, label column):
    /// <ul>
    ///
    /// <li> 1 line per sensor and 3 columns (EEG sensors OR MEG sensors without orientation OR EIT punctual patches)
    ///       <ul TYPE="circle">
    ///       <li> the 1st, 2nd and 3rd columns are respectively position coordinates x, y, z of sensor  </li>
    ///       </ul>
    /// </li>
    /// <li> 1 line per sensor and 4 columns (EEG EIT patches (circular patches)) :
    ///       <ul TYPE="circle">
    ///       <li> the 1st, 2nd and 3rd are respectively position coordinates x, y, z of sensor  </li>
    ///       <li> the 4th is the patche radius (unit relative to the mesh)  </li>
    ///       </ul>
    /// </li>
    /// <li> 1 line per sensor and 6 columns (MEG sensors) :
    ///       <ul TYPE="circle">
    ///       <li> the 1st, 2nd and 3rd are respectively position coordinates x, y, z of sensor  </li>
    ///       <li> the 4th, 5th and 6th are coordinates of vector orientation </li>
    ///       </ul>
    /// </li>
    /// <li> 1 line per integration point for each sensor and 7 columns (MEG sensors) :
    ///       <ul TYPE="circle">
    ///       <li> the 1st, 2nd and 3rd are respectively position coordinates x, y, z of sensor  </li>
    ///       <li> the 4th, 5th and 6th are coordinates of vector orientation </li>
    ///       <li> the 7th is the weight to apply for numerical integration (uses sensor name) </li>
    ///       </ul>
    /// </li>
    /// </ul>

    class OPENMEEG_EXPORT Sensors {

        void init_labels(const Strings& labels) {
            m_pointSensorIdx = std::vector<size_t>(labels.size());
            for (std::size_t i=0; i<labels.size(); ++i)
                m_pointSensorIdx[i] = getSensorIdx(m_names[i]);
        }

    public:

        /// Default constructor. Number of sensors = 0.

        Sensors(): m_nb(0),geometry(nullptr) { }

        /// Default constructor with a geometry. Number of sensors = 0.

        Sensors(const Geometry& g): m_nb(0),geometry(&g) { }

        /// Construct from text file.

        Sensors(const char* filename): geometry(nullptr) { load(filename,'t'); }

        /// Construct from file and geometry (for EIT).

        Sensors(const char* filename,const Geometry& g): geometry(&g) { load(filename,'t'); }

        Sensors(const Matrix& positions,const Geometry& g):
            m_nb(positions.nlin()),m_positions(positions),m_radii(m_nb),geometry(&g)
        {
            m_radii.set(0.0);
            findInjectionTriangles();
        }

        Sensors(const Strings& labels,const Matrix& positions,const Matrix& orientations,const Vector& weights,const Vector& radii):
            m_nb(labels.size()),m_names(labels),m_positions(positions),m_orientations(orientations),m_weights(weights),m_radii(radii)
        {
            init_labels(labels);
        }

        Sensors(const Strings& labels,const Matrix& positions,const Matrix& orientations,const Vector& weights,const Vector& radii,const Geometry& g):
            m_nb(labels.size()),m_names(labels),m_positions(positions),m_orientations(orientations),m_weights(weights),m_radii(radii),geometry(&g)
        {
            // find triangles on which to inject the currents and compute weights
            findInjectionTriangles();
            init_labels(labels);
        }

        /// Load sensors from file. Filetype is 't' for text file or 'b' for binary file.

        void load(const char* filename,const char filetype='t');
        void load(const std::string& filename,const char filetype='t') { load(filename.c_str(),filetype); }

        /// Load description file of sensors from stream. 

        void load(std::istream& in);

        void save(const char* filename) const;
        void save(const std::string& filename) const { save(filename.c_str()); }

        /// Return the number of sensors.

        size_t getNumberOfSensors()   const { return m_nb; }

        /// Return the number of integration points.

        size_t getNumberOfPositions() const { return m_positions.nlin(); }

        /// Return sensors positions.

              Matrix& getPositions()       { return m_positions; }
        const Matrix& getPositions() const { return m_positions; }

        /// Return sensors orientations.

              Matrix& getOrientations()       { return m_orientations; }
        const Matrix& getOrientations() const { return m_orientations; }

        /// Return sensors names.

        Strings& getNames()       { return m_names; }
        Strings  getNames() const { return m_names; }

        bool hasRadii()        const { return m_radii.nlin()>0;        } ///< Return true if contains radii.
        bool hasOrientations() const { return m_orientations.nlin()>0; } ///< Return true if contains orientations.
        bool hasNames()        const { return m_names.size()==m_nb;    } ///< Return true if contains all sensors names.

        /// Return the position (3D point) of the integration point idx.

        Vector getPosition(const size_t idx) const { return m_positions.getlin(idx); }

        /// Return the orientations (3D point) of the integration point idx.

        Vector getOrientation(const size_t idx) const { return m_orientations.getlin(idx); }

        /// Return the name of the idx_th sensor.

        std::string getName(const size_t idx) const {
            om_assert(idx<m_names.size());
            return m_names[idx];
        }

        /// Set the position (3D point) of the integration point i.

        void setPosition(const size_t idx,const Vector& pos) { return m_positions.setlin(idx,pos); }

        /// Set the orientation (3D point) of the integration point i.

        void setOrientation(const size_t idx,const Vector& orient) { return m_orientations.setlin(idx,orient); }

        bool hasSensor(const std::string& name) const { 
            return std::find(m_names.begin(),m_names.end(),name)!=m_names.end();
        }

        size_t getSensorIdx(const std::string& name) const {
            const auto& it = std::find(m_names.begin(),m_names.end(),name);
            if (it==m_names.end()) {
                std::cerr << "Unknown sensor : " << name << std::endl;
                exit(1);
            }
            return std::distance(m_names.begin(),it);
        }

        /// For EIT, get triangles under the current injection electrode.

        Triangles getInjectionTriangles(const size_t idx) const {
            om_assert(idx<m_triangles.size());
            return m_triangles[idx];
        }

        Vector getRadii()   const { return m_radii; }
        Vector getWeights() const { return m_weights; }

        SparseMatrix getWeightsMatrix() const {
            SparseMatrix weight_matrix(getNumberOfSensors(),getNumberOfPositions());
            for(size_t i=0; i<getNumberOfPositions(); ++i)
                weight_matrix(m_pointSensorIdx[i],i) = m_weights(i);
            return weight_matrix;
        }

        /// Return if the sensors object is empty. The sensors object is empty if its number of sensors is null.

        bool isEmpty() { return m_nb==0; }

        /// \brief get info about sensors.

        void info() const;

    private:

        void findInjectionTriangles(); ///< Get the triangles under each EIT sensors.

        size_t                 m_nb;             ///< Number of sensors.
        Strings                m_names;          ///< List of sensors names.
        Matrix                 m_positions;      ///< Matrix of sensors positions. ex: positions(i,j) with  j in {0,1,2} for sensor i.
        Matrix                 m_orientations;   ///< Matrix of sensors orientations. ex: orientation(i,j) with  j in {0,1,2} for sensor i.
        Vector                 m_weights;        ///< Weights of integration points.
        Vector                 m_radii;          ///< Areas of the EIT sensors.
        std::vector<Triangles> m_triangles;      ///< Triangles under each EIT sensors.
        const Geometry*        geometry;         ///< Geometry on which are applied EIT sensors.
        std::vector<size_t>    m_pointSensorIdx; ///< Correspondance between point id and sensor id.
    };
}
