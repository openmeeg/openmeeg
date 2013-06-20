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

#ifndef MESHDESCRIPTION_READER_H
#define MESHDESCRIPTION_READER_H

#include <iostream>
#include <fstream>
#include <string>
#include <stringstream>
#include <vector>
#include <new>
#include <IOUtils.H>

#include <vertex.h>
#include <triangle.h>
#include <meshes.h>
#include <interfaces.h>
#include <domain.h>
#include <MeshDescription/Exceptions.H>

namespace MeshDescription {

    // Parsing of the head description file.
    // See comments for the syntax of the file.
    class Reader {

    public:


    private:

        // Read the mesh section of the description file.
        bool LoadMeshFile(std::istream& is) { //TODO

            bool found;

            is >> io_utils::skip_comments("#") >> io_utils::match_optional("MeshFile", found);

            if ( found ) {

                std::string name;
                is >> io_utils::skip_comments("#") >> io_utils::filename(name, '"', false);
                // Interface::set_file_format(meshes()[i], name);

                //  Load the mesh and check that it is compatible with the first one.
                const std::string& full_name = file_name(name);
                std::ifstream ifs(full_name.c_str());
                if (!ifs.is_open()) {
                    throw MeshDescription::OpenError(full_name);
                }
                load_vtp(ifs);
            }
            return found;
        }

        // Read the interface section of the description file.
        // Check their compatibility and create a data structure indexing all these interfaces.
        void LoadInterfaces(std::istream& is) {
            for (size_t i = 0; i < interfaces().size(); ++i) {
                if (keyword == "Mesh") {
                    std::string name;
                    is >> io_utils::skip_comments("#") >> io_utils::filename(name, '"', false);
                    stringstream ss;
                    ss << i;
                    interfaces()[i].name() = ss.str();
                    meshes[i].mesh_vertices() = vertices; 
                    // Load the mesh
                    const std::string& full_name = file_name(name);
                    std::ifstream ifs(full_name.c_str());
                    if (!ifs.is_open()) {
                        throw MeshDescription::OpenError(full_name);
                    }
                } else if (keyword == "NamedMesh") { // TODO
                    std::string name;
                    is >> io_utils::skip_comments("#") >> io_utils::filename(name, '"', false);
                    stringstream ss;
                    ss << i;
                    interfaces()[i].name() = ss.str();
                    // Load the mesh
                    const std::string& full_name = file_name(name);
                    std::ifstream ifs(full_name.c_str());
                    if (!ifs.is_open()) {
                        throw MeshDescription::OpenError(full_name);
                    }
                } else if (keyword == "Interface") {
                    interfaces()[i].meshes(meshes);
                    std::stringstream ss;
                    ss << i;
                    interfaces()[i].name() = ss.str();
                    getline(ifs, line);
                    std::istringstream iss(line);
                    iss >> interfaces()[i];
                } else { // then NamedInterface
                    interfaces()[i].meshes(meshes);
                    getline(ifs, line);
                    std::istringstream iss(line);
                    std::string i_name;
                    iss >> i_name;
                    if (*i_name.end() == ':') {
                        i_name.substr(0, i_name::npos-1);
                        interfaces()[i].name() = i_name;
                    } else {
                        throw WrongFileFormat(geometry);
                    }
                    iss >> interfaces()[i];
                }
                // unsigned offset_vertices  = vertices.size(); // TODO does it go here ?
                // unsigned offset_triangles = triangles.size();
            }
        }

        // Load the domain part of the description file.
        void LoadDomains(std::istream& is) {

            for (typename Domains::iterator i = domains().begin(); i! = domains().end(); ++i) {

                std::string line;

                is >> io_utils::skip_comments('#') >> io_utils::match("Domain") >> i->name();

                getline(is,line);
                std::istringstream iss(line);
                int number;
                while (iss >> number) {
                    const unsigned index = abs(number);
                    if ((index==0) || (index>interfaces().size())) {
                        throw MeshDescription::NonExistingDomain(i->name(), index);
                    }
                    i->push_back(HalfSpace(number));
                }
                iss.clear();
                iss >> line;
            }
        }

    public:

        Reader(const char* geometry);

              Domains& domains()       { return doms; }
        const Domains& domains() const { return doms; }

        Names domain_names() const {
            Names nms;
            nms.reserve(domains().size());
            for (Domains::const_iterator i = domains().begin(); i!= domains().end(); ++i) {
                nms.push_back(i->name());
            }
            return nms;
        }

        const std::string& domain_name(const unsigned i) const { return domains()[i-1].name(); }

        const Names&       interface_names() const                { return names;      }
        const std::string& interface_name(const unsigned i) const { return names[i-1]; }

              Interfaces& interfaces()       { return interfs; }
        const Interfaces& interfaces() const { return interfs; }

        void stats() const {
             std::cerr << "There are:" << std::endl
                       << interfaces().size() << " interfaces" << std::endl
                       << doms.size()         << " domains"    << std::endl;
        }

    private:

#       if WIN32
        static const char PathSeparator[];
#       else
        static const char PathSeparator   = '/';
#       endif

        static bool IsRelativePath(const std::string& name) {
#           if WIN32
            const char c0 = name[0];
            if (c0=='/' || c0=='\\')
                return false;
            const char c1 = name[1];
            const char c2 = name[2];
            return !(std::isalpha(c0) && c1==':' && (c2=='/' || c2=='\\'));
#           else
            return (name[0]!=PathSeparator);
#           endif
        }

        static std::string find_path(const std::string& name) {
            const std::string::size_type pos = name.find_last_of(PathSeparator);
            return (pos==std::string::npos) ? "" : name.substr(0,pos+1);
        }

        std::string file_name(const std::string& name) const {
            if (!IsRelativePath(name))
                return name;
            return path+name;
        }

        const std::string path;      // Path of the initial geometry file.
        Names             names;     // The file names for the different interfaces.
        Vertices          vertices;  // Vertices descriptions.
        Normals           normals;   // Normals descriptions.
        Triangles         triangles; // Triangles descriptions.
        Meshes            meshes;    // The various meshes composed of triangles.
        Interfaces        interfs;   // The various interfaces between domains composed of meshes.
        Domains           doms;      // Domain descriptions in terms of interfaces.
    };

#   if WIN32
    const char Reader::PathSeparator[] = "/\\";
#   endif

    Reader::Reader(const char* geometry): path(find_path(geometry)) {

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

        std::ifstream ifs(geometry);
        if (!ifs.is_open()) {
            throw MeshDescription::OpenError(geometry);
        }

        int version[2];
        ifs >> io_utils::match("# Domain Description ") >> version[0] >> io_utils::match(".") >> version[1];
        if (ifs.fail()) {
            throw WrongFileFormat(geometry);
        }

        // Process meshes.
        bool meshfile = false;
        if (version[0] == 1) {
            if (version[1] == 0) {
                std::cerr << "Please cpnsider updating the version of the domain description to 1.1 in : " << geometry << std::endl;
            } else if (version[1] == 1) {
                meshfile = LoadMeshFile(ifs);
            } else {
                std::cerr << "Domain Description version not available !" << std::endl;
                throw WrongFileFormat(geometry);
            }
        } else {
            std::cerr << "Domain Description version not available !" << std::endl;
            throw WrongFileFormat(geometry);
        }

        // Process interfaces.
        unsigned num_interfaces;
        std::string InterfaceType;

        ifs >> io_utils::skip_comments('#')
            >> io_utils::match("Interfaces") >> num_interfaces >> InterfaceType;

        if ((InterfaceType == "Mesh")|(InterfaceType == "NamedMesh")|(InterfaceType == "NamedInterface")|(InterfaceType == "Interface")) {
            Interface::keyword = InterfaceType;
        } else {
            throw WrongFileFormat(geometry);
        }

        if (ifs.fail()) {
            throw WrongFileFormat(geometry);
        }

        names.resize(num_interfaces);
        interfaces().resize(num_interfaces);
        if (!meshfile) {
            meshes().resize(num_interfaces);
            size_t i=0;
            for (Interface::iterator mit = interfaces().begin(); mit != interfaces().end(); mit++, i++ ) {
                mit = &meshes()[i];
            }
        }

        LoadInterfaces(ifs);

        // Process domains.
        unsigned num_domains;
        ifs >> io_utils::skip_comments('#') >> io_utils::match("Domains") >> num_domains;

        if (ifs.fail()) {
            throw WrongFileFormat(geometry);
        }

        domains().resize(num_domains);
        LoadDomains(ifs);

        if (ifs.fail()) {
            throw WrongFileFormat(geometry);
        }

        // Close the input file.
        ifs.close();
    }
}

#endif  // ! MESHDESCRIPTION_READER_H
