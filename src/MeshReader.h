/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#ifndef MESHDESCRIPTION_READER_SPECIALIZED_H
#define MESHDESCRIPTION_READER_SPECIALIZED_H

#include <vector>
#include <mesh3.h>
#include <MeshDescription/Reader.H>

namespace MeshReader {

    struct VoidGeometry {
        typedef struct {} Type;
        static void Load(std::istream&) { }
    };

    struct MeshInterface {
        typedef Mesh Type;
        static void set_file_format(Type& interface,const std::string& name) { interface.getFileFormat(name.c_str()); }
        static const char keyword[];
    };
    const char MeshInterface::keyword[] = "Mesh";

    struct Reader: public MeshDescription::Reader<MeshInterface,VoidGeometry> {

        typedef MeshDescription::Reader<MeshInterface,VoidGeometry> base;

    public:

        typedef MeshInterface::Type Mesh;
        typedef base::Domains       Domains;

        Reader(const char* geometry): base(geometry) { }
        
        std::vector<int> sortInterfaceIDAndDomains();
    };

    std::vector<int> Reader::sortInterfaceIDAndDomains() {
#if 0
        for(int i=0;i<doms.size();i++){
            std::cout << "Domaine " << i << " : " << doms[i].name() << std::endl;
            for(int j=0;j<doms[i].size();j++){
                std::cout << "\tInterfaceId : " << doms[i][j].interface() << std::endl;
                std::cout << "\tInOut : " << doms[i][j].inout() << std::endl;
            }
        }
#endif

        // Look for first internal surface :
        int firstSurfId = -1;
        for (Domains::iterator i=domains().begin();i!=domains().end();++i) {
            MeshDescription::Domain& domain = *i;
            if (domain.size()==1 && domain[0].inout()==MeshDescription::Inside) {
                firstSurfId = domain[0].interface();
                std::swap(domain,*domains().begin());
                break;
            }
        }

        assert(firstSurfId >= 0);

        std::vector<int> sortedListOfSurfId;

        // Sort surfaces from first internal to the most external :
        std::vector<bool> domainSeen(domains().size(),false);
        domainSeen[0] = true;
        std::vector<int> sortedDomainsId;
        sortedDomainsId.push_back(0);

        unsigned int currentExternalSurfaceId = firstSurfId;

        for (Interfaces::const_iterator k=interfaces().begin();k!=interfaces().end();++k)
            for (Domains::const_iterator i=domains().begin();i!=domains().end();++i) {
                const unsigned ind = domains().index(*i);
                if (!domainSeen[ind])
                    for (MeshDescription::Domain::const_iterator j=i->begin();j!=i->end();++j) {
                        // find the shared surface which is :
                        //   ** external for the last domain added
                        //   ** internal for the current domain
                        //   ** add the external surface to the list
                        if((j->inout()==MeshDescription::Outside) && (j->interface()==currentExternalSurfaceId)) {
                            sortedListOfSurfId.push_back(currentExternalSurfaceId);
                            sortedDomainsId.push_back(ind);
                            domainSeen[ind] = true;
                            if (i->size()==2)
                                currentExternalSurfaceId = (++j!=i->end()) ? j->interface() : i->begin()->interface();
                        }
                    }
            }

        if (sortedListOfSurfId.size()!=interfaces().size()) {
            std::cout << "Current list : \t" ;
            for(size_t i=0;i<sortedListOfSurfId.size();i++)
                std::cout << sortedListOfSurfId[i] << "\t" ;
            std::cerr << std::endl << "Cannot find " << interfaces().size();
            std::cerr << " nested interfaces with geometry file" << std::endl;
            exit(1);
        }

        // Reordering domains
        std::vector<MeshDescription::Domain> oldDomains = domains();
        for (Domains::iterator i=domains().begin();i!=domains().end();++i)
            *i = oldDomains[sortedDomainsId[domains().index(*i)]];

#if 1
        std::cout << "Sorted List : \t" ;
        for(size_t i=0;i<sortedListOfSurfId.size();i++)
            std::cout << sortedListOfSurfId[i] << " " ;
        std::cout << std::endl;

        std::cout << "Sorted Domains : \t" ;
        for (Domains::const_iterator i=domains().begin();i!=domains().end();++i)
            std::cout << i->name() << "\t";
        std::cout << std::endl;
#endif
        return sortedListOfSurfId;
    }
}


#endif
