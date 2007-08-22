#ifndef MESHDESCRIPTION_READER_SPECIALIZED_H
#define MESHDESCRIPTION_READER_SPECIALIZED_H

#include "mesh3.h"
#include "MeshDescriptionReader.h"
#include "om_utils.h"

namespace MeshDescription {

    struct VoidGeometry {
        typedef struct {} Type;
        static void Load(std::istream&,Type&) { }
    };

    struct MeshInterface {

        typedef Mesh Type;

        struct Validation {
            Validation(const Type& Mesh) { }
            void operator()(const char* name,const Type& interface) { }
        };

        static const char keyword[];

    };
    const char MeshInterface::keyword[] = "Mesh";
//======================================================================================
    template <>
    void Reader<MeshInterface, VoidGeometry>::LoadInterfaces(std::istream& is) {

        // std::cout << "Incoming in MeshdescriptionReaderSpecialized.h - Reader::LoadInterfaces " << std::endl;
        //  The first interface is special as it determines the size of the meshed domain.
        //  Load the first interface and register it for validation.
        std::string interface_name;
        is >> io_utils::skip_comments('#') >> interface_name;

        std::ifstream ifs0(interface_name.c_str());
        if(!ifs0.is_open()) {
            // TODO : Try path relative to position of .geom file to avoid absolute paths in .geom files
            std::cerr << "Error opening file: " << interface_name.c_str() << std::endl;
            exit(1);
        }

        interfaces[0].getFileFormat(interface_name.c_str());
        ifs0 >> interfaces[0];

        MeshInterface::Validation validate(interfaces[0]);
        for (unsigned i=1;i<interfaces.size();++i) {

            is >> io_utils::skip_comments("#") >> interface_name;

            //  Load the interface and check that it is compatible with the first one.
            const char* name = interface_name.c_str();
            std::ifstream ifs(name);
            if(!ifs.is_open()) {
                std::cerr << "Error opening file: " << interface_name.c_str() << std::endl;
                exit(1);
            }

            interfaces[i].getFileFormat(interface_name.c_str());

            ifs >> interfaces[i];
            validate(name,interfaces[i]);
        }
    }


//======================================================================================
    template <>
    std::vector<int> Reader<MeshInterface, VoidGeometry>::sortInterfaceIDAndDomains() {
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
        int innerDomainId;
        for(size_t i=0;i<doms.size();i++)
            if(doms[i].size()==1)
                if(doms[i][0].inout()==0){
                    firstSurfId = doms[i][0].interface();
                    innerDomainId = i;
                    Domain tmpDom;
                    tmpDom = doms[i];
                    doms[i] = doms[0];
                    doms[0] = tmpDom;
                    break;
                }

        assert(firstSurfId >= 0);

        std::vector<int> sortedListOfSurfId;

        // Sort surfaces from first internal to the most external :
        std::vector<bool> domainSeen(doms.size(),false);
        domainSeen[0] = true;
        std::vector<int> sortedDomainsId;
        sortedDomainsId.push_back(0);

        int currentExternalSurfaceId = firstSurfId;

        bool outerSurfaceReached = false;
        for( unsigned int k = 1; k <= interfaces.size(); k += 1 )
        {
            for(size_t i=1;i<doms.size();i++)
            {
                if((domainSeen[i] == false))
                {
                    for(size_t j=0;j<doms[i].size();j++)
                    {
                        // find the shared surface which is :
                        //   ** external for the last domain added
                        //   ** internal for the current domain
                        //   ** add the external surface to the list
                        if((doms[i][j].inout()==1) && (doms[i][j].interface()==currentExternalSurfaceId)){
                            sortedListOfSurfId.push_back(currentExternalSurfaceId);
                            sortedDomainsId.push_back(i);
                            domainSeen[i] = true;

                            if(doms[i].size() == 2){
                                currentExternalSurfaceId = doms[i][(j+1)%2].interface();
                            } else {
                                outerSurfaceReached = true;
                            }
                        }
                    }
                }
            }
        }
        if (sortedListOfSurfId.size() != interfaces.size())
        {
            std::cout << "Current list : \t" ;
            for(size_t i=0;i<sortedListOfSurfId.size();i++)
                std::cout << sortedListOfSurfId[i] << "\t" ;
            std::cerr << std::endl << "Cannot find " << interfaces.size();
            std::cerr << " nested interfaces with geometry file" << std::endl;
            exit(1);
        }

        // Reordering domains
        std::vector<Domain> oldDomains = doms;
        for( unsigned int i = 0; i < sortedDomainsId.size(); i += 1 )
        {
            doms[i] = oldDomains[sortedDomainsId[i]];
        }

#if 0
        std::cout << "Sorted List : \t" ;
        for(int i=0;i<sortedListOfSurfId.size();i++)
            std::cout << sortedListOfSurfId[i] << " " ;
        std::cout << std::endl;

        std::cout << "Sorted Domains : \t" ;
        for(int i=0;i<doms.size();i++)
            std::cout << doms[i].name() << "\t";
        std::cout << std::endl;
#endif
        return sortedListOfSurfId;
    }

}


#endif
