#include <fstream>
#include <cstring>

#include "matrice.h"
#include "IOUtils.h"
#include "sensors.h"

int main(int argc, char** argv){
// usage : sensors sensors_file_description.txt 
	//std::cout << "Number of arguments : " << argc << std::endl;
	//std::cout << "argv[1] = " << argv[1] << std::endl;
	
	int n;
	
#if 0	
	/**** test on empty sensors by default constructor ****/
	sensors emptySensors;
	n = emptySensors.getNumberOfSensors();
	std::cout << "Number of sensors of empty sensors : " << n << std::endl;	
	matrice* p = emptySensors.getSensorsPositions();
	matrice* o = emptySensors.getSensorsOrientations();
	
	if (p==NULL ) std::cout << " p is NULL " << std::endl;
	if (o==NULL ) std::cout << " o is NULL " << std::endl;
	if (emptySensors.isEmpty())
		std::cout << "emptySensors is Empty " << std::endl;
	
	delete p;
	delete o;
#endif	

	
	/*** tests on sensors file ****/
	sensors S(argv[1]);
	
	n = S.getNumberOfSensors();
	std::cout << "Number of sensors of S : " << n << std::endl;
	
	if (S.isEmpty())
		std::cout << "WARNING : empty sensors !" << std::endl;
	else{
		matrice* p = S.getSensorsPositions();
		matrice* o = S.getSensorsOrientations();	
		if(p!=NULL){
			std::cout << std::endl << "Positions of sensors : " << std::endl;
			/*for(int i=0;i<(*p).nlin();i++){
				for(int j=0; j<(*p).ncol();j++)
					std::cout << (*p)(i,j) << "\t" ;
				std::cout << std::endl;
			}*/
		}
		else
			std::cout << "WARNING: no position !" << std::endl;
		if(o!=NULL){
			std::cout << std::endl << "Orientation of sensors : " << std::endl;
			/*for(int i=0;i<(*o).nlin();i++){
				for(int j=0; j<(*o).ncol();j++)
					std::cout << (*o)(i,j) << "\t" ;
				std::cout << std::endl;
			}*/
		}
		else
			std::cout << "WARNING: no orientation !" << std::endl;
	
	
	
		/**** test on copy constructor ****/
		sensors constructorCopyS(S);
		if(constructorCopyS.getNumberOfSensors()!=n)
			std::cout << "ERROR in copy from copy constructor : incorrect number of sensors" << std::endl;
		else
			std::cout << "Number of sensors of copy of S : " << constructorCopyS.getNumberOfSensors() << std::endl;
	
		matrice* pCCS = constructorCopyS.getSensorsPositions();
		matrice* oCCS = constructorCopyS.getSensorsOrientations();	
		
		if(p!=pCCS){ // p and pCCS are not NULL	
			if( ((*p).nlin() != (*pCCS).nlin()) or ((*p).ncol() != (*pCCS).ncol()) )
				std::cout << "ERROR in copy from copy constructor : error in return value of getSensorsPositions()" << std::endl;
			else{
				for(int i=0; i<(*p).nlin();i++)
					for(int j=0;j<(*p).ncol();j++)
						if( (*p)(i,j) != (*pCCS)(i,j) ){
							std::cout << "ERROR in copy from copy constructor : error in return value of getSensorsPositions()" << std::endl;
							break;
						}
			}
		}
		
		if(o!=oCCS){ // o and oCCS are not NULL	
			if( ((*o).nlin() != (*oCCS).nlin()) or ((*o).ncol() != (*oCCS).ncol()) )
				std::cout << "ERROR in copy from copy constructor : error in return value of getSensorsOrientation()" << std::endl;
			else{
				for(int i=0; i<(*o).nlin();i++)
					for(int j=0;j<(*o).ncol();j++)
						if( (*o)(i,j) != (*oCCS)(i,j) ){
							std::cout << "ERROR in copy from copy constructor : error in return value of getSensorsOrientation()" << std::endl;
							break;
						}
			}
		}
	
	}
}
