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
#include <mesh_ios.h>

namespace OpenMEEG {

    /// \brief A class to read geometry and cond file
    class GeometryReader {
    public:

        typedef enum { UNKNOWN_VERSION=-1, VERSION10, VERSION11 } VersionId;

        GeometryReader(Geometry& g): geom(g) {};

        /// \brief read a geometry file
        void read_geom(const std::string&);

        /// \brief read a cond file
        void read_cond(const std::string&);

        VersionId version(const unsigned major,const unsigned minor) const {
            VersionId id = UNKNOWN_VERSION;
            if (major==1) {
                if (minor==0)
                    id = VERSION10;
                if (minor==1)
                    id = VERSION11;
            }
            return id;
        }

        VersionId version() const { return version_id; }

    private:

        VersionId version_id;
        Geometry& geom;
    };

    void GeometryReader::read_geom(const std::string& geometry) {
        // Read the head file description and load the information into the data structures.

        // check file data/README.rst for complete details about the .geom format.
        // The syntax of the head description is a header ("# Domain Description 1.1") followed
        // by two or three sections:
        //
        //     - The first section is optional (for backward compatibility).
        //       Starting with the keyword "MeshFile", it follows the path to the VTK/vtp file containing the meshes.
        //       OR
        //       Starting with the keyword "Meshes", it follows the number of meshes, and the paths to the meshes
        //       (with the keyword "Mesh").
        //
        //     - the second section is introduced by the keyword "Interfaces" followed by a number
        //       (the number of interfaces) and optionally the keyword "Mesh" (for backward compatibility).
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

        if (!ifs.is_open())
            throw OpenMEEG::OpenError(geometry);

        //  Get the version of the geometry file format.

        unsigned major,minor; ///< version of the domain description
        ifs >> io_utils::match("# Domain Description ") >> major >> io_utils::match(".") >> minor;

        if (ifs.fail())
            throw OpenMEEG::WrongFileFormat(geometry);

        version_id = version(major,minor);
        if (version_id==VERSION10)
            std::cerr << "(DEPRECATED) Please consider updating your geometry file to the new format 1.1 (see data/README.rst): "
                      << geometry << std::endl;

        if (version_id==UNKNOWN_VERSION) {
             std::cerr << "Domain Description version not available !" << std::endl;
             throw OpenMEEG::WrongFileFormat(geometry);
        }

        // Extract the absolute path of geometry file

        const std::string& path = absolute_path(geometry);

        // Process meshes.

        if (version_id==VERSION11) {
            // Read the mesh section of the description file.
            // Try to load the meshfile (VTK::vtp file)
            // or try to load the meshes
            bool Is_MeshFile, Is_Meshes;
            ifs >> io_utils::skip_comments("#") >> io_utils::match_optional("MeshFile", Is_MeshFile);
            ifs >> io_utils::skip_comments("#") >> io_utils::match_optional("Meshes", Is_Meshes);
            if (Is_MeshFile) {
                std::string name;
                ifs >> io_utils::skip_comments("#") >> io_utils::filename(name, '"', false);
                const std::string& full_name = (is_relative_path(name)) ? path+name : name;
                geom.load_vtp(full_name);
            } else if (Is_Meshes) {
                unsigned nb_meshes;
                ifs >> nb_meshes;
                Meshes meshes(nb_meshes);
                for (unsigned i=0; i<nb_meshes; ++i) {
                    bool unnamed;
                    ifs >> io_utils::skip_comments("#") >> io_utils::match_optional("Mesh:", unnamed);
                    std::string meshname;
                    std::string filename;
                    if (unnamed) {
                        ifs >> io_utils::filename(filename,'"',false);
                        std::stringstream defaultname;
                        defaultname << i+1;
                        meshname = defaultname.str();
                    } else {
                        ifs >> io_utils::match("Mesh") 
                            >> io_utils::token(meshname, ':') 
                            >> io_utils::filename(filename,'"',false);
                    }
                    const std::string& fullname = (is_relative_path(filename)) ? path+filename : filename;

                    // Load the mesh.

                    meshes[i].load(fullname,false);
                    meshes[i].name() = meshname;
                }

                // Now properly load the meshes into the geometry (no duplicated vertices)

                geom.import_meshes(meshes);
            }
        }

        // Process interfaces.

        unsigned nb_interfaces;
        bool trash; // backward compatibility, catch "Mesh" optionnally.
        ifs >> io_utils::skip_comments('#')
            >> io_utils::match("Interfaces") >> nb_interfaces >> io_utils::match_optional("Mesh", trash);

        if (ifs.fail())
            throw OpenMEEG::WrongFileFormat(geometry);

        // Load interfaces

        std::string id; // id of mesh/interface/domain
        Interfaces interfaces;

        // If meshes are not already loaded

        if (geom.meshes().size()==0) {
            geom.meshes().reserve(nb_interfaces);

            struct MeshDescription {
                std::string interfacename;
                MeshIO*     io;
            };  

            std::vector<MeshDescription> mesh_descriptions(nb_interfaces);

            // First read the total number of vertices

            for (unsigned i=0; i<nb_interfaces; ++i) {
                bool unnamed;
                ifs >> io_utils::skip_comments("#") >> io_utils::match_optional("Interface:", unnamed);
                std::string filename;
                std::string interfacename;
                if (unnamed) {
                    ifs >> io_utils::filename(filename,'"',false);
                    std::stringstream defaultname;
                    defaultname << i+1;
                    interfacename = defaultname.str();
                } else if (version_id==VERSION10) { // backward compatibility
                    std::stringstream defaultname;
                    defaultname << i+1;
                    interfacename = defaultname.str();
                    ifs >> io_utils::filename(filename,'"',false);
                } else {
                    ifs >> io_utils::match("Interface") 
                        >> io_utils::token(interfacename, ':') 
                        >> io_utils::filename(filename,'"',false);
                }
                const std::string& fullname = (is_relative_path(filename)) ? path+filename : filename;
                MeshIO* io = MeshIO::create(fullname);
                io->open();
                io->load_points(geom); 
                mesh_descriptions[i] = MeshDescription({interfacename,io});
            }

            // Second really load the meshes

            for (const auto& desc : mesh_descriptions) {
                geom.meshes().emplace_back(geom,desc.interfacename);
                Mesh& mesh = geom.meshes().back();
                desc.io->load_triangles(mesh);
                interfaces.push_back(Interface(desc.interfacename));
                Interface& interface = interfaces.back();
                interface.oriented_meshes().push_back(OrientedMesh(mesh));
            }
        } else {
            for (unsigned i=0; i<nb_interfaces; ++i) {
                bool unnamed;
                std::string line; // extract a line and parse it
                ifs >> io_utils::skip_comments("#");
                std::getline(ifs, line);
                std::istringstream iss(line);
                iss >> io_utils::match_optional("Interface:", unnamed);
                std::string interfacename;
                if (unnamed) {
                    std::stringstream defaultname;
                    defaultname << i+1;
                    interfacename = defaultname.str();
                } else {
                    iss >> io_utils::match("Interface")
                        >> io_utils::token(interfacename, ':');
                }
                interfaces.push_back(interfacename);
                while (iss >> id) {
                    OrientedMesh::Orientation oriented = OrientedMesh::Normal; // does the id starts with a '-' or a '+' ?
                    if ((id[0]=='-') || (id[0]=='+')) {
                        oriented = (id[0]=='+') ? OrientedMesh::Normal : OrientedMesh::Opposite;
                        id = id.substr(1,id.size());
                    }
                    interfaces[i].oriented_meshes().push_back(OrientedMesh(geom.mesh(id),oriented));
                }
            }
        }

        // Process domains.

        unsigned num_domains;
        ifs >> io_utils::skip_comments('#') >> io_utils::match("Domains") >> num_domains;

        if (ifs.fail())
            throw OpenMEEG::WrongFileFormat(geometry);

        geom.domains().resize(num_domains);
        for (auto& domain : geom.domains()) {
            std::string line;
            ifs >> io_utils::skip_comments('#') >> io_utils::match("Domain");
            if (version_id==VERSION10) { // backward compatibility
                ifs >> domain.name();
            } else {
                ifs >> io_utils::token(domain.name(),':');
            }
            getline(ifs,line);
            std::istringstream iss(line);
            while (iss >> id) {
                bool found = false;
                SimpleDomain::Side side = SimpleDomain::Outside; // does the id starts with a '-' or a '+' ?
                if ((id[0]=='-') || (id[0]=='+')) {
                    side = (id[0]=='-') ? SimpleDomain::Inside : SimpleDomain::Outside;
                    id = id.substr(1,id.size());
                } else if (id=="shared") {
                    std::cerr << "(DEPRECATED) Keyword shared is useless. Please consider updating your geometry file to the new format 1.1 (see data/README.rst): " << geometry << std::endl;
                    break;                    
                }
                for (auto& omesh : interfaces)
                    if (omesh.name()==id) {
                        found = true;
                        if (!omesh.is_mesh_orientations_coherent()) { // check and correct global orientation
                            std::cerr << "Interface \"" << omesh.name() << "\" is not closed !" << std::endl;
                            std::cerr << "Please correct a mesh orientation when defining the interface in the geometry file." << std::endl;
                            exit(1);
                        }
                        domain.boundaries().push_back(SimpleDomain(omesh,side));
                    }

                if (!found)
                    throw OpenMEEG::NonExistingDomain<std::string>(domain.name(), id);
            }
        }

        if (ifs.fail())
            throw OpenMEEG::WrongFileFormat(geometry);
    }

    void GeometryReader::read_cond(const std::string& condFileName) {

        typedef Utils::Properties::Named<std::string,Conductivity<double>> HeadProperties;
        HeadProperties properties(condFileName.c_str());

        // Store the internal conductivity of the external boundary of domain i
        // and store the external conductivity of the internal boundary of domain i

        for (auto& domain : geom.domains()) {
            try {
                const Conductivity<double>& cond = properties.find(domain.name());
                domain.set_conductivity(cond.sigma());
            } catch (const Utils::Properties::UnknownProperty<HeadProperties::Id>& e) {
                throw OpenMEEG::BadDomain(domain.name());
            }
        }
    }
}
