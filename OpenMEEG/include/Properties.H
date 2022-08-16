// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <string>
#include <map>
#include <fstream>
#include <stdexcept>

#include <IOUtils.H>
#include <DataTag.H>

namespace Utils {

    namespace Properties {

        class BadPropertyFile: public std::runtime_error {
        public:
            BadPropertyFile(const std::string& name):
                std::runtime_error("Wrong property file"),filename(name) { }
            ~BadPropertyFile() noexcept { }
            std::string filename;
        };

        template <typename ID>
        class UnknownProperty: public std::runtime_error {
        public:
            UnknownProperty(const ID& i): std::runtime_error("Unknown property"), id(i) { }
            ~UnknownProperty() throw() { }
            ID id;
        };

        template <typename ID>
        class MultiplyDefinedProperty: public std::runtime_error {
        public:
            MultiplyDefinedProperty(const ID& i): std::runtime_error("Property defined multiple times"), id(i) { }
            ~MultiplyDefinedProperty() throw() { }
            ID id;
        };

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

            try {
                Types::DataTag<PROPERTIES> tag;
                is >> io_utils::match("# Properties Description 1.0 (") >> tag >> io_utils::match(")");
                if (!is)
                    throw BadPropertyFile(loader.filename);
            } catch(...) {
                throw BadPropertyFile(loader.filename);
            }
            
            while (is.peek()!=std::istream::traits_type::eof()) {

                Id       id;
                Property property;

                is >> io_utils::skip_comments("#") >> id >> property >> std::ws;
                try {
                    loader.properties.define(id,property);
                } catch(MultiplyDefinedProperty<Id>& e) {
                    std::cerr << "Double definition of " << e.id << " in " << loader.filename
                              << "( error " << e.what() << ")." << std::endl
                              << "Keeping first definition..." << std::endl;
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

            const Property& find(const Id& id) const {
                auto it = base::find(id);
                if (it==base::end())
                    throw UnknownProperty<Id>(id);
                return it->second; 
            }

            void define(const Id& id,const Property& prop) {
                if (base::find(id)!=base::end())
                    throw MultiplyDefinedProperty<Id>(id);
                base::insert(typename base::value_type(id,prop));
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
