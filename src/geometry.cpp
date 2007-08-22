#include "vect3.h"
#include "triangle.h"
#include "mesh3.h"
#include "geometry.h"

#include "MeshDescriptionReaderSpecialized.h"
#include "PropertiesSpecialized.h"


int Geometry::read(char* geomFileName, char* condFileName){

    destroy();

    int result = 0;

    typedef MeshDescription::Reader<MeshDescription::MeshInterface,MeshDescription::VoidGeometry> Reader;
    Reader reader(geomFileName);

    std::vector<Mesh>& Meshes = reader.getInterfaces();
    std::vector<int> meshOrder = reader.sortInterfaceIDAndDomains();

    n = Meshes.size();
    M = new Mesh[n];

    for (int i=0; i<n; i++ ) {
        M[i]=Meshes[meshOrder[i]];
    }

    for (int i=0;i<n;i++) {
        M[i].make_links();
        result += M[i].nbPts();
        result += M[i].nbTrgs();
    }

    std::cout << "Total number of points and triangles : " << result << std::endl;

    std::vector<std::string> domainNames = reader.getDomainNames();

    typedef Utils::Properties::Named< std::string , Conductivity<double> > HeadProperties;
    HeadProperties properties(condFileName);

    sigin = new double[n];
    sigout = new double[n];

    // Store the internal conductivity
    const Conductivity<double>& cond_init=properties.find(domainNames[0]);
    sigin[0] = cond_init.sigma();

    // Store the internal conductivity of the external boundary of domain i
    // and store the external conductivity of the internal boundary of domain i
    for(size_t i=1;i<domainNames.size()-1;i++){
        const Conductivity<double>& cond=properties.find(domainNames[i]);
        sigin[i] = cond.sigma();
        sigout[i-1] = sigin[i];
    }

    const Conductivity<double>& cond_final=properties.find(domainNames[domainNames.size()-1]);
    sigout[n-1] = cond_final.sigma();


    std::cout << "\nChecking" << std::endl;
    for(int i=0;i<n;i++)
        std::cout << "\tMesh " << i << " : internal conductivity = " << sigin[i] << " and external conductivity = " << sigout[i] << std::endl;

    return result;
}


