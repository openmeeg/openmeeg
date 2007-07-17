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
// #include "symmatrice.h" is needed because of warning "files in symmatrice_dcl.h are used but never defined. Indeed, some declaration  in symmatrice_dcl.h don't have definition in matrice.h but in symmatrice.h


// Doc can be generated with doxygen

/*!
 *  Sensors class for EEG and MEG sensors. 
 *  This class is made for read sensors description file. This description file is a file text which can take the shape of :
 *  <ul>
 *	<li> 1 line per 1 sensor and 7 columns (MEG sensors) : 
 *		<ul TYPE="circle">
 *		<li> the 1st column is sensors names </li>
 *		<li> the 2nd, 3rd and 4th are respectively positions coordinates x, y, z of sensor  </li>
 *		<li> the 5th, 6th and 7th are coordinates of vector orientation </li>
 *		</ul>
 *  </li>
 *  <li> 1 line per 1 sensor and 6 columns (MEG sensors) :
 *		<ul TYPE="circle">
 *		<li>- the 1st, 2nd and 3rd are respectively positions coordinates x, y, z of sensor  </li>
 *		<li>- the 4th, 5th and 6th are coordinates of vector orientation </li>
 *		</ul>
 *  </li>
 *	<li> 1 line per 1 sensor and 4 columns (EEG sensors or MEG sensors without orientation) :
 *		<ul TYPE="circle">
 *		<li>- the 1st column is sensors names </li>
 *		<li>- the 2nd, 3rd and 4th are respectively positions coordinates x, y, z of sensor  </li>
 *		</ul>
 *  </li>
 *	<li> 1 line per 1 sensor and 4 columns (EEG sensors or MEG sensors without orientation) :
 *		<ul TYPE="circle">
 *		<li>- the 1st, 2nd and 3rd are respectively positions coordinates x, y, z of sensor  </li>
 *		</ul>
 *  </li>
 *  </ul>
 */
 
 
class sensors{
private:
	int nb_; 							/*!< Number of sensors. */ 
	std::vector<std::string> id_;		/*!< List of sensors name. */ 
	matrice* positions_;				/*!< Pointer on matrix of sensors positions. ex: *positions(i,j) is the jeme coordinate (j in {0,1,2} is {x,y,z}) of the ieme sensor (i=0,..,nb-1). */
	matrice* orientations_; 			/*!< Pointer on matrix of sensors orientations. ex: *orientation(i,j) is the jeme coordinate (j in {0,1,2} is {x,y,z}) of the ieme sensor orientation vector. */

	void copy(const sensors& S); 		/*!< Copy function. Copy sensor S in current sensor object. ex. senors S1; ...; sensors S2(S1); */
	
public:
	sensors(): nb_(0), positions_(NULL), orientations_(NULL) {}  	/*!< Default constructor. Number of sensors = 0. */
	sensors(char* filename, char filetype = 't'); 					/*!< Constructor from path of sensors description file. Option 't' is for text file, and 'b' is for binary file. */
	sensors(const sensors& S) { *this = S; } 						/*!< Copy constructor. */
	~sensors(){
		nb_=0;
		delete positions_;
		delete orientations_;
	} 																/*!< Destructor. Number of sensors = 0. */ 
	
	sensors& operator=(const sensors& S);							 /*!< Copy operator. Copy sensor S in current sensor object. ex. sensors S1; ...; sensors S2 = S1; */
	
	void load(char* filename, char filetype = 't' ); 				/*!< Load description file of sensors from path. Filetype is 't' for text file or 'b' for binary file. */
	void load(std::istream &in); 									/*!< Load description file of sensors from stream. */
	
	int getNumberOfSensors() { return this->nb_; }					/*!< Return the number of sensors. */
	
	matrice* getSensorsPositions() { return this->positions_ ; }		/*!< Return the pointer on matrix of sensors positions. If it is NULL, there are not position. */
	matrice* getSensorsOrientations() {return this->orientations_ ; } 	/*!< Return the pointer on matrix of sensors orientation. If it is NULL, ther are not orientation. */
	std::vector<std::string> getIdOfSensors() {return this->id_ ; }		/*!< Return the vector of whole sensors names. */
		
	int getIndexOfId(std::string id );			/*!< Return the index of id string looked up. */
	vect3 getPosition(std::string id );			/*!< Return the position (3D point) of the id string looked up. */
	vect3 getOrientation(std::string id );		/*!< Return the orientation vector (3D point) of the id string looked up .*/
	
	bool isEmpty() { if(nb_ == 0) return true; else return false; }		/*!< Retunr whether the sensors object is empty. The sensors object is empty if its number of sensors is null. */
};






sensors::sensors(char* filename, char filetype){
	this->load(filename,'t');
}

void sensors::copy(const sensors& S){
	nb_ = S.nb_;
	if ( nb_ != 0 ) {
	    if ( S.id_.size() != 0 ){
	    	for( int i=0; i<nb_; i++)
	    		id_.push_back(S.id_[i]);	
	    }
		
		if ( S.positions_ != NULL ){
			positions_ = new matrice( S.positions_->nlin(), S.positions_->ncol() );
			for( int i=0; i<S.positions_->nlin(); i++)
				for( int j=0; j<S.positions_->ncol(); j++)
					(*positions_)(i,j) = (*S.positions_)(i,j);
		} 
		else positions_ = NULL;
		
		if ( S.orientations_ != NULL ){
			orientations_ = new matrice( S.orientations_->nlin(), S.orientations_->ncol() );
			for( int i=0; i<S.orientations_->nlin(); i++)
				for( int j=0; j<S.orientations_->ncol(); j++)
					(*orientations_)(i,j) = (*S.orientations_)(i,j);
		} 
		else orientations_ = NULL;		
	} 
}

sensors& sensors::operator=(const sensors& S){
	if ( this != &S ) copy(S);
	return *this;
}


void sensors::load(std::istream &in){

	in >> io_utils::skip_comments('#');
    
    std::string s;
    std::getline(in, s);  

	std::stringstream is(s);
	std::string buf; 
    std::vector<std::string> tokens; 

    std::vector<vect3> positions;
    std::vector<vect3> orientations;
    
	// Get data type : 
    int num_of_columns=0;
    while( is >> buf )
    	num_of_columns++;
    
    in.seekg (0, std::ios::beg); // move the get pointer to the beginning of the file.
    
    if (num_of_columns == 7 ){ // id, positions and orientations have to be read.
    	while( !in.eof() ){
    		// get element from line in tokens vector :
    		std::getline(in,s);
    		if(s== "" ) break;
    		std::stringstream iss(s);

    		tokens.clear(); 
       		while( iss >> buf)
        		tokens.push_back(buf);
        	
           	// convert coordonates read as string to double and store them in coord vector :
        	std::vector<double> coord;
    		for(int i=0; i<6; i++){ 
    			std::stringstream tmp_is(tokens[i+1]);
    			double tmp;
    			tmp_is >> tmp;
    			coord.push_back(tmp);
    		}
    		// set position 3D vector :
    		vect3 p(coord[0],coord[1],coord[2]); 
    		positions.push_back(p); 
    		
    		// set orientation 3D vector :
    		vect3 o(coord[3],coord[4],coord[5]);
			orientations.push_back(o);    	
			
			// set private id 
			id_.push_back(tokens[0]);
    	}
    }
    else{
	    if (num_of_columns == 6){ //  positions and orientations have to be read. id not exist.
	    	while( !in.eof()){
	    		// get element from line in tokens vector :
	    		std::getline(in,s);
	    		if(s== "" ) break;
    			std::stringstream iss(s);
    			
    			//tokens.clear();
	    		while( iss >> buf )
	        		tokens.push_back(buf);
	        	
	        	// convert coordonates read as string to double and store them in coord vector :
	        	std::vector<double> coord;
	    		for(int i=0; i<6; i++){ 
	    			std::stringstream tmp_is(tokens[i]);
	    			double tmp;
	    			tmp_is >> tmp;
	    			coord.push_back(tmp);
	    		}
	    		
	    		// set position 3D vector :
	    		vect3 p(coord[0],coord[1],coord[2]); 
	    		positions.push_back(p); 
	    		
	    		// set orientation 3D vector :
	    		vect3 o(coord[3],coord[4],coord[5]);
				orientations.push_back(o);  	
	    	}
	    }
    	else{
			if (num_of_columns == 4){ // id and positions have to be read. orientations not exist.
				while( !in.eof()){
					// get element from line in tokens vector :
					std::getline(in,s);
					if(s== "" ) break;
    				std::stringstream iss(s);
    				
    				//tokens.clear();
					while( iss >> buf )
			    		tokens.push_back(buf);
			    	
			    	id_.push_back(tokens[0]);
			    	
			    	// convert coordonates read as string to double and store them in coord vector :
			    	std::vector<double> coord;
					for(int i=0; i<3; i++){ 
						std::stringstream tmp_is(tokens[i+1]);
						double tmp;
						tmp_is >> tmp;
						coord.push_back(tmp);
					}
					
					// set position 3D vector :
					vect3 p(coord[0],coord[1],coord[2]); 
					positions.push_back(p); 
					
					// set private id 
					id_.push_back(tokens[0]);   	
				}   
			}
    		else{
			    if (num_of_columns == 3){ // only positions have to be read.
			    	while( !in.eof()){
			    		// get element from line in tokens vector :
			    		std::getline(in,s);
			    		if(s== "" ) break;
    					std::stringstream iss(s);
    					
    					//tokens.clear();
			    		while( iss >> buf )
			        		tokens.push_back(buf);

			        	// convert coordonates read as string to double and store them in coord vector :
			        	std::vector<double> coord;
			    		for(int i=0; i<3; i++){ 
			    			std::stringstream tmp_is(tokens[i]);
			    			double tmp;
			    			tmp_is >> tmp;
			    			coord.push_back(tmp);
			    		}
			    		
			    		// set position 3D vector :
			    		vect3 p(coord[0],coord[1],coord[2]); 
			    		positions.push_back(p);  		        	
			    	}   
			    }
			    else{
			    	std::cout << "ERROR: wrong input file" << std::endl;
			    	exit(1);	
			    }
			} // end of else{ if nun_of_colums == 3 }
		} // end of else{if num_of_colums == 4 }
	} // ens of else{if num_of_colums == 6 }


	// init private members :
	nb_ = positions.size();
		
	positions_ = new matrice( nb_, 3);
	for(int i=0; i<positions_->nlin(); i++){
		(*positions_)(i,0) = (positions[i])._x(); 
		(*positions_)(i,1) = (positions[i])._y(); 
		(*positions_)(i,2) = (positions[i])._z(); 
	}
	
	if(!orientations.empty()){
		orientations_ = new matrice( nb_,3);
		for(int i=0; i<positions_->nlin(); i++){
			(*orientations_)(i,0) = (orientations[i])._x();
			(*orientations_)(i,1) = (orientations[i])._y();
			(*orientations_)(i,2) = (orientations[i])._z();
		}
	}	
		
}


void sensors::load(char* filename, char filetype){
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
    sensors::load(in);
    in.close();
}


int sensors::getIndexOfId(std::string id ){
	int i=0;
	while(i<id_.size() && id_[i]!=id)
		i++;
	if(id_[i]!=id)
		return i;
	else
		{ std::cout <<"ERROR: this id not exist! " << std::endl; exit(1); }
}

vect3 sensors::getPosition(std::string id ){
	int ind = getIndexOfId(id);
	vect3 pos( (*positions_)(ind,0), (*positions_)(ind,1), (*positions_)(ind,2) );
	return pos;
}

vect3 sensors::getOrientation(std::string id ){
	int ind = getIndexOfId(id);
	vect3 orient( (*orientations_)(ind,0), (*orientations_)(ind,1), (*orientations_)(ind,2) );
	return orient;
}

#endif
