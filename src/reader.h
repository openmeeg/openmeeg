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

#ifndef OPENMEEG_READER_H
#define OPENMEEG_READER_H

#include <vector>
#include <interface.h>
#include <MeshDescription/Exceptions.H>
#include <IOUtils.H>
#include "PropertiesSpecialized.h"

namespace OpenMEEG {

    bool Geometry::is_relative_path(const std::string& name) {
    #if WIN32
        const char c0 = name[0];
        if (c0 == '/' || c0 == '\\') {
            return false;
        }
        const char c1 = name[1];
        const char c2 = name[2];
        return !(std::isalpha(c0) && c1==':' && (c2=='/' || c2=='\\'));
    #else
        return ( name[0] != PathSeparator );
    #endif
    }

    void Geometry::read_geom(const std::string& geometry) {
        // Read the head file description and load the information into temporary data structures.

        // The syntax of the head description is a header ("# Domain Description (1.0):") followed
        // by three sections:
        //
        //     - The first section is made of two fields and defines the geometry of the rectangular
        //       domain. First the origin of the domain is given by the keyword "Origin" and the
        //       vector specifying the coordinates of the upper corner of the domain. Then, the size
        //       of the domain is given by another vector and introduced by the keyword "DomainSize".
        //
        //     - the second section is introduced by the keyword "Interfaces" followed by a number
        //       (the number of interfaces) and a type (currently Mesh, NamedMesh, Interface, NamedInterface are possible).
        //       The section just contains a list of names (one per line, the remainder of the line
        //       being ignored).
        //
        //     - the third section is introduced by the keyword "Domains" and the number of domains
        //       (everything else on the line containing the keyword is currently ignored). The section
        //       contains domains descriptions, one per line. Each domain consist of:
        //
        //         o a domain name.
        //         o a list of signed numbers: the absolute value of the number gives describes an
        //           interface by its index in the "Interfaces" list (indices are starting at one so
        //           that the sign is meaningful), the sign of the number depicts whether the interior
        //           or the exterior of the interface should be used to select the domain.
        //
        // Any line starting with # is considered a comment and is silently ignored.
        
        std::ifstream ifs(geometry.c_str());

        if (!ifs.is_open()) {
            throw MeshDescription::OpenError(geometry);
        }

        unsigned version[2];
        ifs >> io_utils::match("# Domain Description ") >> version[0] >> io_utils::match(".") >> version[1];
        if (ifs.fail()) {
            throw MeshDescription::WrongFileFormat(geometry);
        }

        // extract the absolut path of geometry file
        const std::string::size_type pos = geometry.find_last_of(this->PathSeparator);
        const std::string path = (pos == std::string::npos) ? "" : geometry.substr(0, pos+1);

        // Process meshes. -----------------------------------------------------------------------------------
        bool meshfile = false;
        if (version[0] == 1) {
            if (version[1] == 0) {
                std::cerr << "Please consider updating the version of the domain description to 1.1 in the geometry file: "
                          << geometry << std::endl;
            } else if (version[1] == 1) {
                // Read the mesh section of the description file.
                // Try to load the meshfile (VTK::vtp file)
                bool meshfile;
                ifs >> io_utils::skip_comments("#") >> io_utils::match_optional("MeshFile", meshfile);
                if ( meshfile ) {
                    std::string name;
                    ifs >> io_utils::skip_comments("#") >> io_utils::filename(name, '"', false);
                    //  Load the mesh and check that it is compatible with the first one.
                    const std::string& full_name = (is_relative_path(name))?path+name:name;
                    std::ifstream ifs(full_name.c_str());
                    if (!ifs.is_open()) {
                        throw MeshDescription::OpenError(full_name);
                    }
                    // load_vtp_file(ifs);
                }
            } else {
                std::cerr << "Domain Description version not available !" << std::endl;
                throw MeshDescription::WrongFileFormat(geometry);
            }
        } else {
            std::cerr << "Domain Description version not available !" << std::endl;
            throw MeshDescription::WrongFileFormat(geometry);
        }

        // Process interfaces. -----------------------------------------------------------------------------------
        unsigned num_interfaces;
        std::string interfaceType;

        ifs >> io_utils::skip_comments('#')
            >> io_utils::match("Interfaces") >> num_interfaces >> interfaceType;

        if (ifs.fail()) {
            throw MeshDescription::WrongFileFormat(geometry);
        }

        // load the interfaces/meshes
        std::string id; // id of mesh/interface/domain
        Interfaces interf;
        if ( interfaceType == "Mesh"||"NamedMesh") { // ---------------------------------------
            meshes().reserve(num_interfaces);
            std::string interfacename[num_interfaces], filename[num_interfaces], fullname[num_interfaces]; // names
            // First read the total number of vertices
            unsigned nb_vertices = 0;
            for (unsigned i = 0; i < num_interfaces; i++ ) {
                if (interfaceType == "Mesh") {
                    ifs >> io_utils::skip_comments("#") >> io_utils::filename(filename[i], '"', false);
                    std::stringstream defaultname;
                    defaultname << i+1;
                    interfacename[i] = defaultname.str();
                } else {
                    ifs >> io_utils::skip_comments("#") >> interfacename[i] >> io_utils::filename(filename[i], '"', false);
                    if (*interfacename[i].end() == ':') {
                        interfacename[i] = interfacename[i].substr(0, interfacename[i].npos-1);
                    } else {
                        throw MeshDescription::WrongFileFormat(geometry);
                    }
                }
                Mesh m;
                fullname[i] = (is_relative_path(filename[i]))?path+filename[i]:filename[i];
                nb_vertices += m.load(fullname[i].c_str(), false, false); 
            }
            vertices().reserve(nb_vertices);
            // Second load the mesh
            for (unsigned i = 0; i < num_interfaces; i++ ) {
                Mesh m(vertices(), interfacename[i]);
                meshes().push_back(m);
                meshes()[i].load(fullname[i].c_str());
                interf.push_back( Interface(interfacename[i]) );
                interf[i].push_back(&(meshes()[i])); // one mesh per interface: mesh at this adress
            }
        } else if (interfaceType == "Interface"||"NamedInterface") { // -----------------------
            std::string interfacename;
            for (unsigned i = 0; i < num_interfaces; i++ ) {
                std::string line; // extract a line and parse it
                std::getline(ifs, line);
                std::istringstream iss(line);
                if (interfaceType == "Interface") {
                    std::stringstream defaultname;
                    defaultname << i+1;
                    interfacename = defaultname.str();
                } else {
                    iss >> interfacename;
                    if (*interfacename.end() == ':') {
                        interfacename = interfacename.substr(0, interfacename.npos-1);
                    } else {
                        throw MeshDescription::WrongFileFormat(geometry);
                    }
                }
                interf.push_back( interfacename );
                while (iss >> id) {
                    interf[i].push_back(&mesh(id));
                }
            }
        } else {
            throw MeshDescription::WrongFileFormat(geometry);
        }

        // Process domains. -----------------------------------------------------------------------------------
        unsigned num_domains;
        ifs >> io_utils::skip_comments('#') >> io_utils::match("Domains") >> num_domains;

        if (ifs.fail()) {
            throw MeshDescription::WrongFileFormat(geometry);
        }

        domains().resize(num_domains);
        for (Domains::iterator dit = domain_begin(); dit != domain_end(); ++dit) {
            std::string line;
            ifs >> io_utils::skip_comments('#') >> io_utils::match("Domain") >> dit->name();
            getline(ifs, line);
            std::istringstream iss(line);
            while (iss >> id) {
                bool found = false;
                for (Interfaces::iterator iit = interf.begin(); iit != interf.end() ; iit++) {
                    bool inside = (id[0] == '-'); // does the id starts with a '-' ?
                    if (iit->name() == (inside?id.substr(1, id.npos):id)) { // TODO (+)
                        found = true;
                        dit->push_back(HalfSpace(*iit, inside));
                    }
                }
                if (!found) {
                    throw MeshDescription::NonExistingDomain(dit->name(), 0); //TODO I don't want to give 0 index but name ! template Exceptions?
                }
            }
        }

        // Search for the innermost (resp. outermost) domain and set a boolean on the domain in the vector domains.
        // An innermost (resp. outermost) domain is (here) defined as the only domain (inside/outside) represented by only one interface.
        // Specials domain names can enforce which domain is innermost ("Innermost") or outermost: ("Outermost", "Air")
        Domains::iterator dit_out;
        bool outer;
        for (Domains::iterator dit = domain_begin(); dit != domain_end(); ++dit) {
            outer = true;
            for (Domain::iterator hit = dit->begin(); hit != dit->end(); ++hit) {
                outer = outer && !(hit->inside());
                dit_out = dit;
            }
            if (outer) {
                dit_out->outermost() = true;
                for (Domain::iterator hit = dit_out->begin(); hit != dit_out->end(); ++hit) {
                    hit->interface().set_to_outermost();
                }
                break;
            }
        }

        if (ifs.fail()) {
            throw MeshDescription::WrongFileFormat(geometry);
        }

        // Close the input file. -----------------------------------------------------------------------------------
        ifs.close();
    }

    void Geometry::read_cond(const std::string condFileName) {

        typedef Utils::Properties::Named< std::string , Conductivity<double> > HeadProperties;
        HeadProperties properties(condFileName.c_str());

        // Store the internal conductivity of the external boundary of domain i
        // and store the external conductivity of the internal boundary of domain i
        for (Domains::iterator dit = this->domain_begin(); dit != domain_end(); dit++) {
            const Conductivity<double>& cond = properties.find(dit->name());
            dit->sigma() =  cond.sigma();
        }
    }
}

#endif  //! OPENMEEG_READER_H
