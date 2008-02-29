/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
lastrevision      : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#ifndef UTILS_PROPERTIES
#define UTILS_PROPERTIES

#include <string>
#include <map>
#include <fstream>

#include <IOUtils.H>
#include <DataTag.H>

namespace Utils {

    namespace Properties {

        template <typename PROPERTIES>
        struct PropertyLoader {

            PropertyLoader(PROPERTIES& prop): filename("[unkmown]"),properties(prop) { }

            PropertyLoader(const char* name,PROPERTIES& prop): filename(name),properties(prop) {
                std::ifstream ifs(filename);
                ifs >> *this;
            }

            const char* filename;
            PROPERTIES& properties;
        };

        template <typename PROPERTIES>
        std::istream& operator>>(std::istream& is,PropertyLoader<PROPERTIES>& loader) {

            typedef typename PROPERTIES::Id       Id;
            typedef typename PROPERTIES::Property Property;

            Types::DataTag<PROPERTIES> tag;

            try {
                is >> io_utils::match("# Properties Description 1.0 (") >> tag >> io_utils::match(")");
            } catch(...) {
                std::cerr << "Wrong property file: " << loader.filename << '.' << std::endl
                          << "  Check whether the file exists and contains the proper properties." << std::endl;
                exit(1);
            }

            while (is.peek()!=std::istream::traits_type::eof()) {

                Id       id;
                Property property;

                is >> io_utils::skip_comments("#") >> id >> property >> std::ws;
                try {
                    loader.properties.define(id,property);
                } catch(int& i) {
                    std::cerr << "Double definition of " << id << " in " << loader.filename << "( error " << i << ")." << std::endl
                              << "Keeping first definition...." << std::endl;
                }
            }
            return is;
        }

        //  A class that can associate various properties to identities.
        //  Can be used as a basis class for converting names into effective values.
        //  This is mainly to be used in the case where there is a small number of
        //  properties.

        template <typename ID,typename PROPERTY>
        class Named: public std::map<ID,PROPERTY> {
            
            typedef std::map<ID,PROPERTY> base;

        public:

            typedef ID       Id;
            typedef PROPERTY Property;

            Named(): base() { }

            Named(const char* filename): base() { PropertyLoader<Named> loader(filename,*this); }

            //  Define and retrieve a property using its id.

            const Property& find(const Id& id) const { return base::find(id)->second; }

            void define(const Id& id,const Property& prop) throw(char*) {
                if (base::find(id)!=base::end()) {
                    throw "Double definition";
                } else {
                    base::insert(typename base::value_type(id,prop));
                }
            }
        };

        template <typename ID,typename PROPERTY>
        std::istream& operator>>(std::istream& is,Named<ID,PROPERTY>& properties) {
            PropertyLoader<Named<ID,PROPERTY> > loader(properties);
            return is >> loader;
        }

        template <typename ID,typename PROPERTY>
        std::ostream& operator<<(std::ostream& os,const Named<ID,PROPERTY>& properties) {
            for (typename Named<ID,PROPERTY>::const_iterator i=properties.begin();i!=properties.end();++i)
                os << i->first << ' ' << i->second << std::endl;
            return os;
        }
    }
}

#endif  //  ! UTILS_PROPERTIES
