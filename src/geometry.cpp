#include "vect3.h"
#include "triangle.h"
#include "mesh3.h"
#include "om_utils.h"
#include "geometry.h"
#include "MeshReader.h"
#include "PropertiesSpecialized.h"

int Geometry::read(const char* geomFileName, const char* condFileName) {

    destroy();

    has_cond = false; // default param

    size_t npts = 0;
    size_t ntrgs = 0;

    MeshReader::Reader reader(geomFileName);

    std::vector<Mesh>& Meshes = reader.interfaces();
    std::vector<int> meshOrder = reader.sortInterfaceIDAndDomains();

    n = Meshes.size();
    M = new Mesh[n];

    for (int i=0;i<n;i++)
        M[i] = Meshes[meshOrder[i]];

    for (int i=0;i<n;i++) {
        M[i].make_links();
        npts += M[i].nbPts();
        ntrgs += M[i].nbTrgs();
    }

    std::cout << "Total number of points    : " << npts << std::endl;
    std::cout << "Total number of triangles : " << ntrgs << std::endl;

    std::vector<std::string> domainNames = reader.domain_names();

    if(condFileName)
    {
        has_cond = true;
        typedef Utils::Properties::Named< std::string , Conductivity<double> > HeadProperties;
        HeadProperties properties(condFileName);

        sigin  = new double[n];
        sigout = new double[n];

        // Store the internal conductivity
        const Conductivity<double>& cond_init=properties.find(domainNames[0]);
        sigin[0] = cond_init.sigma();

        // Store the internal conductivity of the external boundary of domain i
        // and store the external conductivity of the internal boundary of domain i
        for(size_t i=1;i<domainNames.size()-1;i++) {
            const Conductivity<double>& cond=properties.find(domainNames[i]);
            sigin[i] = cond.sigma();
            sigout[i-1] = sigin[i];
        }

        const Conductivity<double>& cond_final=properties.find(domainNames[domainNames.size()-1]);
        sigout[n-1] = cond_final.sigma();

        std::cout << "\nChecking" << std::endl;
        for(int i=0;i<n;i++)
            std::cout << "\tMesh " << i << " : internal conductivity = " << sigin[i] << " and external conductivity = " << sigout[i] << std::endl;
    }

    m_size = npts + ntrgs;
    return m_size;
}

bool Geometry::selfCheck() const {
    bool OK = true;
    for(size_t i = 0; i < nb(); ++i)
    {
        const Mesh& m1 = getM(i);
        if(m1.selfIntersection())
        {
            warning(std::string("Mesh is self intersecting !"));
            m1.info();
            OK = false;
        }
        for(size_t j = i+1; j < nb(); ++j)
        {
            const Mesh& m2 = getM(j);
            if(m1.intersection(m2))
            {
                warning(std::string("2 meshes are intersecting !"));
                m1.info();
                m2.info();
                OK = false;
            }
        }
    }
    return OK;
}

bool Geometry::check(const Mesh& m) const {
    bool OK = true;
    if(m.selfIntersection())
    {
        warning(std::string("Mesh is self intersecting !"));
        m.info();
        OK = false;
    }
    for(size_t i = 0; i < nb(); ++i)
    {
        const Mesh& m1 = getM(i);
        if(m1.intersection(m))
        {
            warning(std::string("Mesh is intersecting with one of the mesh in geom file !"));
            m1.info();
            OK = false;
        }
    }
    return OK;
}
