#ifndef H_sensors
#define H_sensors

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "IOUtils.h"
#include "vect3.h"
#include "matrice.h"
#include "symmatrice.h"

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
 *    <li> 1 line per 1 sensor and 4 columns (EEG sensors or MEG sensors without orientation) :
 *        <ul TYPE="circle">
 *        <li>- the 1st, 2nd and 3rd are respectively positions coordinates x, y, z of sensor  </li>
 *        </ul>
 *  </li>
 *  </ul>
 */

class Sensors {
private:
    size_t m_nb;                    /*!< Number of sensors. */
    std::vector<std::string> m_id;  /*!< List of sensors name. */
    matrice m_positions;            /*!< Matrix of sensors positions. ex: positions(i,j) with  j in {0,1,2} for sensor i */
    matrice m_orientations;         /*!< Matrix of sensors orientations. ex: orientation(i,j) with  j in {0,1,2} for sensor i */

    void copy(const Sensors& S);    /*!< Copy function. Copy sensor S in current sensor object. ex. senors S1; ...; sensors S2(S1); */

public:
    Sensors(): m_nb(0) {} /*!< Default constructor. Number of sensors = 0. */
    Sensors(char* filename, char filetype = 't'); /*!< Construct from file. Option 't' is for text file, and 'b' is for binary file. */
    Sensors(const Sensors& S) { copy(S); }        /*!< Copy constructor. */
    ~Sensors() { m_nb=0; }                        /*!< Destructor. Number of sensors = 0. */

    Sensors& operator=(const Sensors& S); /*!< Copy operator. Copy sensor S in current sensor object. ex. sensors S1; ...; sensors S2 = S1; */

    void load(char* filename, char filetype = 't' ); /*!< Load sensors from file. Filetype is 't' for text file or 'b' for binary file. */
    void load(std::istream &in); /*!< Load description file of sensors from stream. */

    int getNumberOfSensors() const { return m_nb; } /*!< Return the number of sensors. */

    matrice& getSensorsPositions() { return m_positions ; } /*!< Return a reference on sensors positions. */
    matrice getSensorsPositions() const { return m_positions ; } /*!< Return a copy of sensors positions */

    matrice& getSensorsOrientations() {return m_orientations ; } /*!< Return a reference on sensors orientations. */
    matrice getSensorsOrientations() const {return m_orientations ; } /*!< Return a copy of sensors orientations. */

    std::vector<std::string>& getSensorsIds() {return m_id ; } /*!< Return a reference on sensors ids. */
    std::vector<std::string> getSensorsIds() const {return m_id ; } /*!< Return a copy of sensors ids. */

    std::vector<std::string> getIdOfSensors() {return m_id ; } /*!< Return the vector of whole sensors names. */
    bool hasOrientations() const { return m_orientations.nlin() > 0 ;} /*!< Return true if contains orientations */
    bool hasIds() const { return m_id.size() == m_nb ;} /*!< Return true if contains all sensors ids (names) */
    int getIndexOfId(std::string id ); /*!< Return the index of id string looked up. */
    Vect3 getPosition(std::string id ); /*!< Return the position (3D point) of the id string looked up. */
    Vect3 getOrientation(std::string id ); /*!< Return the orientation vector (3D point) of the id string looked up .*/

    bool isEmpty() { if(m_nb == 0) return true; else return false; } /*!< Return if the sensors object is empty. The sensors object is empty if its number of sensors is null. */
};

Sensors::Sensors(char* filename, char filetype) {
    this->load(filename,'t');
}

void Sensors::copy(const Sensors& S) {
    m_nb = S.m_nb;
    if ( m_nb != 0 ) {
        if ( S.m_id.size() != 0 ) {
            for( size_t i=0; i<m_nb; i++)
                m_id.push_back(S.m_id[i]);
        }

        m_positions = matrice( S.m_positions );
        m_orientations = matrice( S.m_orientations );
    }
}

Sensors& Sensors::operator=(const Sensors& S) {
    if ( this != &S ) copy(S);
    return *this;
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
    while( is >> buf )
        num_of_columns++;

    // Determine the number of lines
    in.seekg(0,std::ios::beg);
    size_t num_of_lines = 0;
    while(!in.fail())
    {
        std::getline(in,s);
        num_of_lines++;
    }
    num_of_lines--;

    // init private members :
    m_nb = num_of_lines;
    m_positions = matrice( m_nb, 3);
    if ( num_of_columns > 4 ) {
        m_orientations = matrice( m_nb, 3);
    }

    size_t current_line_id = 0;
    in.clear();
    in.seekg(0,std::ios::beg); // move the get pointer to the beginning of the file.
    while ( !in.fail() ) {
        // Tokenize line
        std::getline(in,s);

        if( s == "" ) break; // Skip blank line
        std::stringstream iss(s);
        tokens.clear();
        while( iss >> buf)
            tokens.push_back(buf);

        assert(tokens.size() == num_of_columns); // Each line has same length

        tokensIterator =  tokens.begin();

        if ( (num_of_columns == 7) or (num_of_columns == 4) ) { // Get label
            m_id.push_back(*tokensIterator);
            tokensIterator++;
        }

        // read position
        std::vector<double> coord;
        for(size_t i=0; i<3; i++){
            std::stringstream tmp_is(*tokensIterator);
            double tmp;
            tmp_is >> tmp;
            m_positions(current_line_id,i) = tmp;
            tokensIterator++;
        }
        
        if ( num_of_columns > 4 ) {
            // read orientation
            std::vector<double> coord;
            for(size_t i=0; i<3; i++){
                std::stringstream tmp_is(*tokensIterator);
                double tmp;
                tmp_is >> tmp;
                m_orientations(current_line_id,i) = tmp;
                tokensIterator++;
            }
        }
        current_line_id++;
    }

}

void Sensors::load(char* filename, char filetype) {
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

int Sensors::getIndexOfId(std::string id ) {
    size_t i=0;
    while(i<m_id.size() && m_id[i]!=id)
        i++;
    if(m_id[i]!=id)
        return i;
    else
        { std::cout <<"ERROR: this id not exist! " << std::endl; exit(1); }
}

Vect3 Sensors::getPosition(std::string id ) {
    int ind = getIndexOfId(id);
    Vect3 pos( m_positions(ind,0), m_positions(ind,1), m_positions(ind,2) );
    return pos;
}

Vect3 Sensors::getOrientation(std::string id ) {
    int ind = getIndexOfId(id);
    Vect3 orient( m_orientations(ind,0), m_orientations(ind,1), m_orientations(ind,2) );
    return orient;
}

inline std::ostream& operator<<(std::ostream& f,const Sensors &S) {
    f << "Nb of sensors : " << S.getNumberOfSensors() << std::endl;
    f << "Positions" << std::endl;
    f << S.getSensorsPositions();
    if(S.hasOrientations())
    {
        f << "Orientations" << std::endl;
        f << S.getSensorsOrientations();
    }
    if(S.hasIds())
    {
        f << "Ids" << std::endl;
        std::vector<std::string> ids = S.getSensorsIds();
        for(size_t i = 0; i < ids.size(); ++i)
        {
            f << ids[i] << std::endl;
        }
    }
    return f;
}

#endif
