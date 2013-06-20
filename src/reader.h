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

namespace OpenMEEG {

    bool Geometry::is_relative_path(const std::string& name) {
        #if WIN32
        const char c0 = name[0];
        if (c0=='/' || c0=='\\')
            return false;
        const char c1 = name[1];
        const char c2 = name[2];
        return !(std::isalpha(c0) && c1==':' && (c2=='/' || c2=='\\'));
        #else
        return (name[0]!=PathSeparator);
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

        int version[2];
        ifs >> io_utils::match("# Domain Description ") >> version[0] >> io_utils::match(".") >> version[1];
        if (ifs.fail()) {
            throw MeshDescription::WrongFileFormat(geometry);
        }

        // extract the absolut path of geometry file
        const std::string::size_type pos = geometry.find_last_of(this->PathSeparator);
        const std::string path = (pos == std::string::npos) ? "" : geometry.substr(0, pos+1);

        // Process meshes.
        bool meshfile = false;
        if (version[0] == 1) {
            if (version[1] == 0) {
                std::cerr << "Please consider updating the version of the domain description to 1.1 in the geometry file: "
                          << geometry << std::endl;
            } else if (version[1] == 1) {
                // Read the mesh section of the description file.
                // Try to load the meshfile ((VTK) vtp file)
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
                    // load_vtp(ifs);
        }
            } else {
                std::cerr << "Domain Description version not available !" << std::endl;
                throw MeshDescription::WrongFileFormat(geometry);
            }
        } else {
            std::cerr << "Domain Description version not available !" << std::endl;
            throw MeshDescription::WrongFileFormat(geometry);
        }

        // Process interfaces.
        unsigned num_interfaces;
        std::string InterfaceType;

        ifs >> io_utils::skip_comments('#')
            >> io_utils::match("Interfaces") >> num_interfaces >> InterfaceType;

        if ((InterfaceType == "Mesh")|(InterfaceType == "NamedMesh")|(InterfaceType == "NamedInterface")|(InterfaceType == "Interface")) {
            Interface::keyword = InterfaceType;
        } else {
            throw MeshDescription::WrongFileFormat(geometry);
        }
        if (ifs.fail()) {
            throw MeshDescription::WrongFileFormat(geometry);
        }

        interfaces().resize(num_interfaces);
        if (!meshfile) {
            meshes().resize(num_interfaces);
            size_t i=0;
            for (Interfaces::iterator iit = interface_begin(); iit != interface_end(); iit++, i++ ) {
                iit->push_back(&(meshes()[i]));
            }
        }

        //  load the interfaces
        for (size_t i = 0; i < interfaces().size(); ++i) {
            if (Interface::keyword == "Mesh") {
                std::string filename;
                ifs >> io_utils::skip_comments("#") >> io_utils::filename(filename, '"', false);
                std::stringstream ss;
                ss << i+1;
                interfaces()[i].name() = ss.str();
                meshes()[i].name()     = ss.str();
                this->vertices().reserve(100);
                // meshes()[i].all_vertices() = &this->vertices();
                // Load the mesh
                const std::string full_name = (is_relative_path(filename))?path+filename:filename;
                mesh_load(full_name.c_str(), meshes()[i]);
            } else if (Interface::keyword == "NamedMesh") { // TODO
                std::string meshname;
                std::string filename;
                ifs >> io_utils::skip_comments("#") >> meshname >> io_utils::filename(filename, '"', false);
                if (*meshname.end() == ':') {
                    meshname = meshname.substr(0, meshname.npos-1);
                } else {
                    throw MeshDescription::WrongFileFormat(geometry);
                }
                interfaces()[i].name() = meshname;
                meshes()[i].name()     = meshname;
                meshes()[i].all_vertices(&vertices()); 
                // Load the mesh
                const std::string full_name = (is_relative_path(filename))?path+filename:filename;
                mesh_load(full_name.c_str(), meshes()[i]);
            } else if (Interface::keyword == "Interface") {
                // interfaces()[i].meshes(meshes);
                std::stringstream ss;
                ss << i+1;
                interfaces()[i].name() = ss.str();
                std::string line;
                std::getline(ifs, line);
                std::istringstream iss(line);
                iss >> interfaces()[i];
            } else { // then NamedInterface
                // interfaces()[i].meshes()=meshes;
                std::string line;
                getline(ifs, line);
                std::istringstream iss(line);
                std::string i_name;
                iss >> i_name;
                if (*i_name.end() == ':') {
                    i_name = i_name.substr(0, i_name.npos-1);
                    interfaces()[i].name() = i_name;
                } else {
                    throw MeshDescription::WrongFileFormat(geometry);
                }
                iss >> interfaces()[i];
            }
        }

        // Process domains.
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
            std::string id; // id = oriented interface name (starting with '-' if denotes inside )
            while (iss >> id) {
                bool found = false;
                for (Interfaces::const_iterator iit = interfaces().begin(); iit != interfaces().end() ; iit++) {
                    bool inside = (id[0] == '-'); // does the id starts with a '-' ?
                    if (iit->name() == (inside?id.substr(1, id.npos):id)) {
                        found = true;
                        dit->push_back(HalfSpace(&*iit, inside));
                    }
                }
                if (!found) {
                    throw MeshDescription::NonExistingDomain(dit->name(), 0); //TODO I don't want to give 0 index but name ! template Exceptions?
                }
            }
            iss.clear();
            iss >> line;
        }

        // Search and place the Air domain as the last domain in the vector // TODO look for the outermost instead of "Air"
        Domains::iterator dit;
        for (dit = domain_begin(); dit != domain_end(); ++dit) {
            if (dit->name() == "Air" ) {
                break;
            }
        }
        std::swap(*dit, domains()[domains().size()-1]);
        for (Domains::const_iterator hit = domain_end(); hit != domain_end(); ++hit) {
            if (domain_end()->size() == 1) {
                for (Domain::const_iterator mit = hit->begin(); mit != hit->end(); mit++ ) {
                    mit->interface()[0]->outermost() = true;
                }
            }
            else {
                throw MeshDescription::NonExistingDomain(dit->name(), 0); //TODO throw or treat, Air domain is defined by several interfaces
            }

        }


        // Search for an innermost domain and place it as the first domain in the vector
        // An innermost domain is (here) defined as the only domain represented by only one interface
        bool only_one = false;
    /*    for (Domains::iterator dit2 = domain_begin(); dit2 != domain_end(); ++dit2) {
            if ( (dit2->size() == 1) && dit2->inside() ) {
                if (only_one) {
                    only_one = false;
                    break;
                }
                else {
                    dit = dit2;
                    only_one = true;
                }
            }
        }
        if (only_one) {
            dit->innermost() = true;
            std::cout << "Innermost domain \"" << dit->name() << "\" found." << std::endl;
            std::swap(*dit, domains()[0]);
        }
*/
        if (ifs.fail()) {
            throw MeshDescription::WrongFileFormat(geometry);
        }

        // Close the input file.
        ifs.close();

        geom_generate_indices();
    }


}

#endif  //! OPENMEEG_READER_H
