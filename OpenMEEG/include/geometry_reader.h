/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

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

#pragma once

#include <vector>
#include <interface.h>
#include <GeometryExceptions.H>
#include <IOUtils.H>
#include <PropertiesSpecialized.h>
#include <geometry_io.h>

namespace OpenMEEG {

    /// \brief A class to read geometry and cond file
    class GeometryReader {
        public:
            GeometryReader(Geometry& g): geom(g) {};

            /// \brief read a geometry file
            void read_geom(const std::string&);

            /// \brief read a cond file
            void read_cond(const std::string&);

        private:
            Geometry& geom;

            /// \return true if name is a realtive path. \param name
            bool is_relative_path(const std::string& name);
        #if WIN32
            static const char PathSeparator[];
        #else
            static const char PathSeparator   = '/';
        #endif
    };
    #if WIN32
        const char GeometryReader::PathSeparator[] = "/\\";
    #endif

    bool GeometryReader::is_relative_path(const std::string& name) {
    #if WIN32
        const char c0 = name[0];
        if ( c0 == '/' || c0 == '\\' ) {
            return false;
        }
        const char c1 = name[1];
        const char c2 = name[2];
        return !( std::isalpha(c0) && c1 == ':' && ( c2 == '/' || c2 == '\\' ) );
    #else
        return ( name[0] != PathSeparator );
    #endif
    }

    void GeometryReader::read_geom(const std::string& geometry) {
        // Read the head file description and load the information into the data structures.

        // check file data/README.rst for complete details about the .geom format.
        // The syntax of the head description is a header ("# Domain Description 1.1") followed
        // by two or three sections:
        //
        //     - The first section is optional (for backward compatibility).
        //       Starting with the keyword "MeshFile", it follows the path to the VTK/vtp file containing the meshes.
        //       OR
        //       Starting with the keyword "Meshes", it follows the number of meshes, and the paths to the meshes (with the keyword "Mesh").
        //
        //     - the second section is introduced by the keyword "Interfaces" followed by a number
        //       (the number of interfaces) and optionnally the keyword "Mesh" (for backward compatibility).
        //       The section contains the list of interfaces preceded by keyword "Interface".
        //
        //     - the third section is introduced by the keyword "Domains" and the number of domains
        //       (everything else on the line containing the keyword is currently ignored). The section
        //       contains domains descriptions, one per line. Each domain consist of:
        //
        //         o a domain name.
        //         o a list of IDs (signed numbers or signed names): the sign ('+'by default) of the ID depicts 
        //           whether the interior or the exterior of the interface should be used to select the domain.
        //
        // Any line starting with # is considered a comment and is silently ignored.
        
        std::ifstream ifs(geometry.c_str());

        if ( !ifs.is_open() ) {
            throw OpenMEEG::OpenError(geometry);
        }

        //  Get the version of the geometry file format.

        unsigned version[2]; ///< version of the domain description
        ifs >> io_utils::match("# Domain Description ") >> version[0] >> io_utils::match(".") >> version[1];

        if ( ifs.fail() ) {
            throw OpenMEEG::WrongFileFormat(geometry);
        }

        geom.version_id = Geometry::UNKNOWN_VERSION;
        if (version[0]==1) {
            if (version[1]==0) {
                geom.version_id = Geometry::VERSION10;
                std::cerr << "(DEPRECATED) Please consider updating your geometry file to the new format 1.1 (see data/README.rst): "
                          << geometry << std::endl;
            }
            if (version[1]==1)
                geom.version_id = Geometry::VERSION11;
        }

        if (geom.version_id==Geometry::UNKNOWN_VERSION) {
             std::cerr << "Domain Description version not available !" << std::endl;
             throw OpenMEEG::WrongFileFormat(geometry);
        }

        // extract the absolut path of geometry file
        const std::string::size_type pos = geometry.find_last_of(this->PathSeparator);
        const std::string path = (pos == std::string::npos) ? "" : geometry.substr(0, pos+1);

        // Process meshes. -----------------------------------------------------------------------------------
        if (geom.version_id==Geometry::VERSION11) {
            // Read the mesh section of the description file.
            // Try to load the meshfile (VTK::vtp file)
            // or try to load the meshes
            bool Is_MeshFile, Is_Meshes;
            ifs >> io_utils::skip_comments("#") >> io_utils::match_optional("MeshFile", Is_MeshFile);
            ifs >> io_utils::skip_comments("#") >> io_utils::match_optional("Meshes", Is_Meshes);
            if ( Is_MeshFile ) {
                std::string name;
                ifs >> io_utils::skip_comments("#") >> io_utils::filename(name, '"', false);
                const std::string& full_name = (is_relative_path(name))?path+name:name;
                geom.load_vtp(full_name);
            } else if ( Is_Meshes ) {
                unsigned nb_meshes;
                ifs >> nb_meshes;
                std::vector<std::string> meshname(nb_meshes); // names
                std::vector<std::string> filename(nb_meshes);
                std::vector<std::string> fullname(nb_meshes);
                Meshes meshes(nb_meshes);
                for ( unsigned i = 0; i < nb_meshes; ++i ) {
                    bool unnamed;
                    ifs >> io_utils::skip_comments("#") >> io_utils::match_optional("Mesh:", unnamed);
                    if ( unnamed ) {
                        ifs >> io_utils::filename(filename[i], '"', false);
                        std::stringstream defaultname;
                        defaultname << i+1;
                        meshname[i] = defaultname.str();
                    } else {
                        ifs >> io_utils::match("Mesh") 
                            >> io_utils::token(meshname[i], ':') 
                            >> io_utils::filename(filename[i], '"', false);
                    }
                    fullname[i] = (is_relative_path(filename[i]))?path+filename[i]:filename[i];
                    // Load the mesh.
                    meshes[i].load(fullname[i], false);
                    meshes[i].name() = meshname[i];
                }
                // Now properly load the meshes into the geometry (not dupplicated vertices)
                geom.import_meshes(meshes);
            }
        }

        // Process interfaces. -----------------------------------------------------------------------------------
        unsigned nb_interfaces;
        bool trash; // backward compatibility, catch "Mesh" optionnally.
        ifs >> io_utils::skip_comments('#')
            >> io_utils::match("Interfaces") >> nb_interfaces >> io_utils::match_optional("Mesh", trash);

        if ( ifs.fail() ) {
            throw OpenMEEG::WrongFileFormat(geometry);
        }

        // load the interfaces
        std::string id; // id of mesh/interface/domain
        Interfaces interfaces;
        // if meshes are not already loaded
        if ( geom.nb_meshes() == 0 ) { // ---------------------------------------
            geom.meshes_.reserve(nb_interfaces);
            std::vector<std::string> interfacename(nb_interfaces);
            std::vector<std::string> filename(nb_interfaces);
            std::vector<std::string> fullname(nb_interfaces);

            // First read the total number of vertices

            unsigned nb_vertices = 0;
            for ( unsigned i = 0; i < nb_interfaces; ++i ) {
                bool unnamed;
                ifs >> io_utils::skip_comments("#") >> io_utils::match_optional("Interface:", unnamed);
                if ( unnamed ) {
                    ifs >> io_utils::filename(filename[i], '"', false);
                    std::stringstream defaultname;
                    defaultname << i+1;
                    interfacename[i] = defaultname.str();
                } else if (geom.version_id==Geometry::VERSION10) { // backward compatibility
                    std::stringstream defaultname;
                    defaultname << i+1;
                    interfacename[i] = defaultname.str();
                    ifs >> io_utils::filename(filename[i], '"', false);
                } else {
                    ifs >> io_utils::match("Interface") 
                        >> io_utils::token(interfacename[i], ':') 
                        >> io_utils::filename(filename[i], '"', false);
                }
                Mesh m;
                fullname[i] = (is_relative_path(filename[i]))?path+filename[i]:filename[i];
                nb_vertices += m.load(fullname[i], false, false); 
            }
            geom.vertices_.reserve(nb_vertices);
            // Second really load the meshes
            for ( unsigned i = 0; i < nb_interfaces; ++i ) {
                geom.meshes_.push_back(Mesh(geom.vertices_, interfacename[i]));
                geom.meshes_[i].load(fullname[i], false);
                interfaces.push_back( Interface(interfacename[i]) );
                interfaces[i].push_back(OrientedMesh(geom.meshes_[i], true)); // one mesh per interface, (well oriented)
            }
        } else { // -----------------------
            std::string interfacename;
            for ( unsigned i = 0; i < nb_interfaces; ++i ) {
                bool unnamed;
                std::string line; // extract a line and parse it
                ifs >> io_utils::skip_comments("#");
                std::getline(ifs, line);
                std::istringstream iss(line);
                iss >> io_utils::match_optional("Interface:", unnamed);
                if ( unnamed ) {
                    std::stringstream defaultname;
                    defaultname << i+1;
                    interfacename = defaultname.str();
                } else {
                    iss >> io_utils::match("Interface")
                        >> io_utils::token(interfacename, ':');
                }
                interfaces.push_back( interfacename );
                while ( iss >> id ) {
                    bool oriented = true; // does the id starts with a '-' or a '+' ?
                    if ( ( id[0] == '-' ) || ( id[0] == '+' ) ) {
                        oriented = ( id[0] == '+' );
                        id = id.substr(1, id.size());
                    }
                    interfaces[i].push_back(OrientedMesh(geom.mesh(id), oriented));
                }
            }
        }

        // Process domains. -----------------------------------------------------------------------------------

        unsigned num_domains;
        ifs >> io_utils::skip_comments('#') >> io_utils::match("Domains") >> num_domains;

        if ( ifs.fail() ) {
            throw OpenMEEG::WrongFileFormat(geometry);
        }

        geom.domains_.resize(num_domains);
        for ( Domains::iterator dit = geom.domain_begin(); dit != geom.domain_end(); ++dit) {
            std::string line;
            if (geom.version_id==Geometry::VERSION10) { // backward compatibility
                ifs >> io_utils::skip_comments('#') >> io_utils::match("Domain") >> dit->name();
            } else {
                ifs >> io_utils::skip_comments('#') >> io_utils::match("Domain") >> io_utils::token(dit->name(), ':');
            }
            getline(ifs, line);
            std::istringstream iss(line);
            while ( iss >> id ) {
                bool found = false;
                bool inside = false; // does the id starts with a '-' or a '+' ?
                if ( ( id[0] == '-' ) || ( id[0] == '+' ) ) {
                    inside = ( id[0] == '-' );
                    id = id.substr(1, id.size());
                } else if ( id == "shared" ) {
                    std::cerr << "(DEPRECATED) Keyword shared is useless. Please consider updating your geometry file to the new format 1.1 (see data/README.rst): " << geometry << std::endl;
                    break;                    
                }
                for ( Interfaces::iterator iit = interfaces.begin(); iit != interfaces.end() ; ++iit) {
                    if ( iit->name() == id ) {
                        found = true;
                        if ( !iit->check() ) { // check and correct global orientation
                            std::cerr << "Interface \"" << iit->name() << "\" is not closed !" << std::endl;
                            std::cerr << "Please correct a mesh orientation when defining the interface in the geometry file." << std::endl;
                            exit(1);
                        }
                        dit->push_back(HalfSpace(*iit, inside));
                    }
                }
                if ( !found ) {
                    throw OpenMEEG::NonExistingDomain<std::string>(dit->name(), id);
                }
            }
        }

        // Search for the outermost domain and set boolean OUTERMOST on the domain in the vector domains.
        // An outermost domain is (here) defined as the only domain outside represented by only one interface.

        Domains::iterator dit_out;
        for ( Domains::iterator dit = geom.domain_begin(); dit != geom.domain_end(); ++dit) {
            bool outer = true;
            for ( Domain::iterator hit = dit->begin(); hit != dit->end(); ++hit)
                outer = outer && !(hit->inside());

            if (outer) {
                dit_out = dit;
                dit_out->outermost() = true;
                for (Domain::iterator hit = dit_out->begin(); hit != dit_out->end(); ++hit) {
                    hit->interface().set_to_outermost();
                }
                break;
            }
        }

        // Determine if the geometry is nested or not
        // The geometry is considered non nested if (at least) one domain is defined as being outside two or more interfaces
        // OR
        // if 2 interfaces are composed by a same mesh oriented once correctly once wrongly.

        bool nested = true;
        for ( Domains::const_iterator dit = geom.domain_begin(); dit != geom.domain_end(); ++dit) {
            unsigned out_interface = 0;
            if ( dit != dit_out ) {
                for ( Domain::const_iterator hit = dit->begin(); hit != dit->end(); ++hit) {
                    if ( !hit->inside() ) {
                        out_interface++;
                    }
                }
            }
            if ( out_interface >= 2 ) {
                nested = false;
                break;
            }
        }

        if ( nested ) {
            for ( Geometry::const_iterator mit = geom.begin(); mit != geom.end(); ++mit) {
                unsigned m_oriented = 0;
                for ( Domains::const_iterator dit = geom.domain_begin(); dit != geom.domain_end(); ++dit) {
                    for ( Domain::const_iterator hit = dit->begin(); hit != dit->end(); ++hit) {
                        for ( Interface::const_iterator iit = hit->first.begin(); iit != hit->first.end(); ++iit) {
                            if ( iit->mesh() == *mit ) {
                                m_oriented += (iit->orientation());
                            }
                        }
                    }
                }
                if ( m_oriented == 0 ) {
                    nested = false; // TODO unless a mesh is defined but unused ...
                    break;
                }
            }
        }
        geom.is_nested_ = nested;

        if ( ifs.fail() ) {
            throw OpenMEEG::WrongFileFormat(geometry);
        }

        // Close the input file. -----------------------------------------------------------------------------------
        ifs.close();
    }

    void GeometryReader::read_cond(const std::string& condFileName) {

        typedef Utils::Properties::Named<std::string, Conductivity<double> > HeadProperties;
        HeadProperties properties(condFileName.c_str());

        // Store the internal conductivity of the external boundary of domain i
        // and store the external conductivity of the internal boundary of domain i
        for ( Domains::iterator dit = geom.domain_begin(); dit != geom.domain_end(); ++dit) {
            try {
                const Conductivity<double>& cond = properties.find(dit->name());
                dit->sigma() =  cond.sigma();
            } catch( const Utils::Properties::UnknownProperty<HeadProperties::Id>& e) {
                throw OpenMEEG::BadDomain(dit->name());
            }
        }
    }
}
