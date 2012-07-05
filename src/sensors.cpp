/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
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

#include "sensors.h"

namespace OpenMEEG {

    void Sensors::copy(const Sensors& S) {
        m_nb = S.m_nb;
        if ( m_nb != 0 ) {
            if ( S.m_names.size() != 0 ) {
                for( size_t i=0; i<m_nb; i++)
                    m_names.push_back(S.m_names[i]);
            }

            m_positions = Matrix( S.m_positions );
            m_orientations = Matrix( S.m_orientations );
        }
    }

    Sensors::Sensors(const char* filename) {
        this->load(filename,'t');
    }

    Sensors& Sensors::operator=(const Sensors& S) {
        if ( this != &S ) copy(S);
        return *this;
    }

    bool Sensors::hasSensor(std::string name) {
        for(size_t i = 0; i < m_names.size(); ++i) {
            if(m_names[i] == name) {
                return true;
            }
        }
        return false;
    }

    size_t Sensors::getSensorIdx(std::string name) {
        for(size_t i = 0; i < m_names.size(); ++i) {
            if(m_names[i] == name) {
                return i;
            }
        }
        std::cerr << "Unknown sensor : " << name << std::endl;
        exit(1);
        return 0;
    }

    void Sensors::load(const char* filename, char filetype) {
        std::ifstream in;
        if(filetype == 't')
            in.open(filename,std::ios::in);
        else
            if(filetype == 'b')
                in.open(filename,std::ios::in|std::ios::binary);
            else
                { std::cout << "ERROR: unkown filetype. " << std::endl; exit(1); }


        if(!in.is_open())
            { std::cerr<<"Error Reading File : " << filename << std::endl;  exit(1);  }
        Sensors::load(in);
        in.close();
    }

    void Sensors::load(std::istream &in) {

        in >> io_utils::skip_comments('#');

        std::string buf;
        std::vector<std::string> tokens;
        std::vector<std::string>::const_iterator tokensIterator;

        // Get data type :
        std::string s;
        std::getline(in, s);
        std::stringstream is(s);
        size_t num_of_columns = 0;
        while(is >> buf)
            num_of_columns++;

        // Determine the number of lines
        in.seekg(0,std::ios::beg);
        size_t num_of_lines = 0;
        while (std::getline(in,s))
            if (!s.empty())
                ++num_of_lines;

        // init private members :
        m_positions = Matrix(num_of_lines,3);
        m_weights = Vector(num_of_lines);
        m_pointSensorIdx = std::vector<size_t>(num_of_lines);

        if (num_of_columns>4)
            m_orientations = Matrix(num_of_lines,3);

        m_nb = 0;
        size_t current_line_id = 0;
        in.clear();
        in.seekg(0,std::ios::beg); // move the get pointer to the beginning of the file.
        while (std::getline(in,s))
            if (!s.empty()) {

                //  Tokenize the line.

                std::stringstream iss(s);
                tokens.clear();
                while (iss >> buf)
                    tokens.push_back(buf);

                if(tokens.size()!= num_of_columns) {
                    std::cout << tokens.size() << " != " << num_of_columns << std::endl;
                    std::cerr << "Problem while reading Sensors file" << std::endl;
                    std::cerr << "Each line should have the same number of elements" << std::endl;
                    exit(1);
                }

                tokensIterator =  tokens.begin();

                size_t sensor_idx = m_nb;
                // See if it's actually a new sensor or just a new integration point
                if ((num_of_columns >= 7) || (num_of_columns == 4)) { // Get label
                    std::string sensor_name = *tokensIterator;
                    tokensIterator++;
                    if(hasSensor(sensor_name)) {
                        sensor_idx = getSensorIdx(sensor_name);
                    } else {
                        m_nb++;
                        m_names.push_back(sensor_name);
                    }
                } else {
                    m_nb++;
                }

                m_pointSensorIdx[current_line_id] = sensor_idx;

                // read position
                for(size_t i=0;i<3;++i){
                    std::stringstream tmp_is(*tokensIterator);
                    double tmp;
                    tmp_is >> tmp;
                    m_positions(current_line_id,i) = tmp;
                    tokensIterator++;
                }

                if (num_of_columns>4) {
                    // read orientation
                    std::vector<double> coord;
                    for(size_t i=0;i<3;++i){
                        std::stringstream tmp_is(*tokensIterator);
                        double tmp;
                        tmp_is >> tmp;
                        m_orientations(current_line_id,i) = tmp;
                        tokensIterator++;
                    }
                }

                // Try to read weight
                if (tokensIterator!=tokens.end()) {
                    std::stringstream ss(*tokensIterator);
                    tokensIterator++;
                    double tmp;
                    ss >> tmp;
                    m_weights(current_line_id) = tmp;
                } else {
                    m_weights(current_line_id) = 1.0;
                }

                assert(tokensIterator==tokens.end()); // Check if everything has been read

                current_line_id++;
            }
    }

    void Sensors::save(const char* filename) {
        std::ofstream outfile(filename);
        for(size_t i = 0; i < getNumberOfPositions(); ++i) {
            if (hasNames())
                outfile << m_names[m_pointSensorIdx[i]] << " ";
            outfile << m_positions.getlin(i) << " ";
            if (hasOrientations())
                outfile << m_orientations.getlin(i) << " ";
            outfile << m_weights(i) << std::endl;
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
}
