#ifndef MESHDESCRIPTION_READER_H
#define MESHDESCRIPTION_READER_H

//#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <new>
#include "IOUtils.h"
#include "mesh3.h"

namespace MeshDescription {

    typedef enum { Inside, Outside } InOut;

    //  Parsing of the head description file.
    //  See comments for the syntax of the file.

    template <typename INTERFACE,typename GEOMETRY>
    class Reader {
        
        //  A domain is the association of a name and a vector of pairs.
        //  The first element of each pair corresponds to an interface.
        //  The second element of each pair states whether the domain is contains in the
        //  inside or the ouside of the volume bounded by the interface.

        typedef unsigned InterfaceId;

        struct ElementaryDomain: protected std::pair<InterfaceId,InOut> {
            typedef std::pair<InterfaceId,InOut> base;
            ElementaryDomain(const int num): base(abs(num)-1,((num<0) ? Inside : Outside)) { }
            InterfaceId interface() const { return base::first;  }
            InOut       inout()     const { return base::second; }
        };

        struct Domain: public std::vector<ElementaryDomain> {
            Domain() { }

            //  The name of the domain.

                  std::string& name()       { return name_; }
            const std::string& name() const { return name_; }

        private:

            std::string name_;       // Name of the domain.
        };

        //  Domains is just a collection of Domain (here a simple vector).

        struct Domains: public std::vector<Domain> {
            unsigned index(const Domain& dom) const { return &dom-&*this->begin(); }
        };

        typedef std::vector<typename INTERFACE::Type> Interfaces;

        typedef typename GEOMETRY::Type Geometry;

        //  Read the interface section of the description file.
        //  Check their compatibility and create a data structure indexing all these interfaces.

        void LoadInterfaces(std::istream& is) {
            
            //  The first interface is special as it determines the size of the meshed domain.
            //  Load the first interface and register it for validation.

            std::string interface_name;
            is >> io_utils::skip_comments('#') >> interface_name;
            std::ifstream ifs0(interface_name.c_str());
            if(!ifs0.is_open()) {
                std::cerr << "Error opening file: " << interface_name.c_str() << std::endl;
                exit(1);
            }
            ifs0 >> interfaces[0];

            typename INTERFACE::Validation validate(interfaces[0]);
            for (unsigned i=1;i<interfaces.size();++i) {

                is >> io_utils::skip_comments("#") >> interface_name;

                //  Load the interface and check that it is compatible with the first one.
                const char* name = interface_name.c_str();
                std::ifstream ifs(name);
                if(!ifs.is_open()) {
                    std::cerr << "Error opening file: " << interface_name.c_str() << std::endl;
                    exit(1);
                }

                ifs >> interfaces[i];
                validate(name,interfaces[i]);
            }
        }

        //  Load the domain part of the description file.

        void LoadDomains(std::istream& is) {

            for (typename Domains::iterator i=domains().begin();i!=domains().end();++i) {

                std::string line;

                is >> io_utils::skip_comments('#') >> io_utils::match("Domain") >> i->name();

                getline(is,line);
                std::istringstream iss(line);
                int number;
                while (iss >> number) {
                    const unsigned index = abs(number);
                    if ((index==0) || (index>interfaces.size()))
                        throw 2;
                    i->push_back(ElementaryDomain(number));
                }
                iss.clear();
                iss >> line;
            }
        }

    public:

        Reader(const char* geometry);
        Reader(const Reader& reader);
        Reader& operator=(Reader& reader);
        
              Domains& domains()       { return doms; }
        const Domains& domains() const { return doms; }

        void stats() const {
            std::cerr << "There are:" << std::endl
                      << interfaces.size() << " interfaces" << std::endl
                      << doms.size()       << " domains"    << std::endl;
        }
        

        
        Interfaces& getInterfaces() { return this->interfaces; }
        std::vector<int> sortInterfaceIDAndDomains(); // implemented in MeshDescriptionReaderSpecialized.h
        std::vector<std::string> getDomainNames() {
            std::vector<std::string> domainNames;
            for(int i=0;i<doms.size();i++)
                domainNames.push_back(doms[i].name());
            return domainNames;
        };
        
    private :    
        Interfaces interfaces;  //  The various levelsets depicting interfaces between domains.
        Domains    doms;        //  Domain descriptions in terms of interfaces.
        Geometry   geom;        //  The geometry of the depicted domain (origin and size).
    };

    template <typename INTERFACES,typename GEOMETRY>
    Reader<INTERFACES,GEOMETRY>::Reader(const char* geometry) {
        //  Read the head file description and load the information into temporary data structures.

        //  The syntax of the head description is a header ("# Domain Description (1.0):") followed
        //  by three sections:
        //
        //      - The first section is made of two fields and defines the geometry of the rectangular
        //        domain. First the origin of the domain is given by the keyword "Origin" and the
        //        vector specifying the coordinates of the upper corner of the domain. Then, the size
        //        of the domain is given by another vector and introduced by the keyword "DomainSize".
        //
        //      - the second section is introduced by the keyword "Interfaces" followed by a number
        //        (the number of interfaces) and a type (currently only "LevelSet" is possible).
        //        The section just contains a list of names (one per line, the remainder of the line
        //        being ignored).
        //
        //      - the third section is introduced by the keyword "Domains" and the number of domains
        //        (everything else on the line containing the keyword is currently ignored). The section
        //        contains domains descriptions, one per line. Each domain consist of:
        //
        //          o a domain name.
        //          o a list of signed numbers: the absolute value of the number gives describes an
        //            interface by its index in the "Interfaces" list (indices are starting at one so
        //            that the sign is meaningful), the sign of the number depicts whether the interior
        //            or the exterior of the interface should be used to select the domain.
        //
        //  Any line starting with # is considered a comment and is silently ignored.
        
        std::ifstream ifs(geometry);
        ifs >> io_utils::match("# Domain Description 1.0");
        GEOMETRY::Load(ifs,geom);
        

        //  Process interfaces.

        unsigned num_interfaces;
        ifs >> io_utils::skip_comments('#')
            >> io_utils::match("Interfaces") >> num_interfaces >> io_utils::match(INTERFACES::keyword);
        if (ifs.fail())
            throw "Wrong file format!";

        interfaces.reserve(num_interfaces);
        interfaces.resize(num_interfaces);  
        LoadInterfaces(ifs);

        //  Process domains.
        
        unsigned num_domains;
        ifs >> io_utils::skip_comments('#') >> io_utils::match("Domains") >> num_domains;
        if (ifs.fail())
            throw "Wrong file format!";
        
        doms.reserve(num_domains);
        doms.resize(num_domains);
        LoadDomains(ifs);
        
        if (ifs.fail())
            throw "Wrong file format!";

        //  Close the input file.

        ifs.close();
    }
#if 0   
    // constructeur de copie :
    template <typename INTERFACES,typename GEOMETRY>
    Reader<INTERFACES,GEOMETRY>::Reader(const Reader<INTERFACES,GEOMETRY>& reader) {
        *this = reader;
    }
    
    // operateur de copie :
    template <typename INTERFACES,typename GEOMETRY>
    Reader<INTERFACES,GEOMETRY>& Reader<INTERFACES,GEOMETRY>::operator=(Reader<INTERFACES,GEOMETRY>& reader){
        if(this == &reader) return reader;
        interfaces = reader.interfaces;
        doms = reader.doms;
        geom = reader.geom;
        return *this;
    }
#endif
}

#endif  // ! MESHDESCRIPTION_READER_H
