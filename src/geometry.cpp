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

#include "vect3.h"
#include "triangle.h"
#include "mesh3.h"
#include "om_utils.h"
#include "geometry.h"
#include "MeshReader.h"
#include "PropertiesSpecialized.h"

namespace OpenMEEG {

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
        for(int i = 0; i < nb(); ++i)
        {
            const Mesh& m1 = getM(i);
            if(m1.selfIntersection())
            {
                warning(std::string("Mesh is self intersecting !"));
                m1.info();
                OK = false;
		std::cout << "Self intersection for mesh number " << i << std:: endl;
            }
            for(int j = i+1; j < nb(); ++j)
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
        for(int i = 0; i < nb(); ++i)
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

    int Geometry::getDomain(const Vect3& p) {
        for (int i=0;i<nb();i++)
            if ((this->getM(i)).containsPoint(p))
                return i;
            else if (i==nb()-1)
                return nb();
        return -1; // should never be reached
    }

}
