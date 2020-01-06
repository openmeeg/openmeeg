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

#include <map>
#include <vector>
#include <string>

#include <interface.h>
#include <GeometryExceptions.H>
#include <IOUtils.H>
#include <GeometryIO.h>
#include <MeshIO.h>

namespace OpenMEEG::GeometryIOs {

    // Read the head file description (.geom files) and load the information into a Geometry.

    // Check file data/README.rst for complete details about the .geom format.
    // The syntax of the head description is a header ("# Domain Description 1.1") followed
    // by two or three sections:
    //
    //     - The first section is optional (for backward compatibility).
    //       Starting with the keyword "MeshFile", it follows the path to the VTK/vtp file containing the meshes.
    //       OR
    //       Starting with the keyword "Meshes", it follows the number of meshes, and the paths to the meshes
    //       (with the keyword "Mesh").
    //
    //     - The second section is introduced by the keyword "Interfaces" followed by a number
    //       (the number of interfaces) and optionally the keyword "Mesh" (for backward compatibility).
    //       The section contains the list of interfaces preceded by keyword "Interface".
    //
    //     - The third section is introduced by the keyword "Domains" and the number of domains
    //       (everything else on the line containing the keyword is currently ignored). The section
    //       contains domains descriptions, one per line. Each domain consist of:
    //
    //         o a domain name.
    //         o a list of IDs (signed numbers or signed names): the sign ('+'by default) of the ID depicts 
    //           whether the interior or the exterior of the interface should be used to select the domain.
    //
    // Any line starting with # is considered a comment and is silently ignored.

    /// \brief Reader for geometry file (.geom).

    class OPENMEEG_EXPORT GeomFile: public GeometryIO {

        typedef GeometryIO base;

        void load_meshes(Geometry& geometry) override;

        typedef enum { UNKNOWN_VERSION=-1, VERSION10, VERSION11 } VersionId;

        const char* name() const override { return "geom"; }
        Matrix load_data() const override { return Matrix(); }

        virtual void save_geom(const Geometry& geometry) override;
        virtual void save_data(const Geometry&,const Matrix&) const override { }
        virtual void write() const override { }

        GeometryIO* clone(const std::string& filename) const override { return new GeomFile(filename); }

        void load_domains(Geometry& geometry) override;

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

        static std::string default_name(const unsigned n) {
            std::stringstream res;
            res << n;
            return res.str();
        }

        std::string given_name(const std::string& keyword) {
            std::string res;
            fs >> io_utils::match(keyword) >> io_utils::token(res,':');
            return res;
        }

        std::string filename() {
            std::string filename;
            fs >> io_utils::filename(filename,'"',false);
            return (is_relative_path(filename)) ? directory+filename : filename;
        }

        static bool extract_sign(std::string& str) {
            const char sign = str[0];
            if (sign=='-' || sign=='+') {
                str = str.substr(1,str.size());
                return sign=='+';
            }
            return true;
        }

        std::string section_name(const unsigned n,const std::string& keyword) {
            bool unnamed;
            fs >> io_utils::skip_comments("#") >> io_utils::match_optional(keyword+':',unnamed);
            const std::string& name = (unnamed || version_id==VERSION10) ? default_name(n+1) : given_name(keyword);
            return name;
        }

        void read_mesh_descriptions(const unsigned n,const std::string& keyword) {
            for (unsigned i=0; i<n; ++i) {
                const std::string& name = section_name(i,keyword);
                const std::string& path = filename();
                meshes.push_back({ name, path });
            }
        }

        //  Extract the tokens on the remainder of the line in the stream.

        std::vector<std::string> get_line_tokens() {
            std::string line;
            std::getline(fs,line);
            std::istringstream iss(line);
            typedef std::istream_iterator<std::string> Iterator;
            std::vector<std::string> res;
            std::copy(Iterator(iss),Iterator(),std::back_inserter(res));
            return res;
        }

        GeomFile(const std::string& filename=""): base(filename,"geom") { }

        static const GeomFile prototype;

        typedef Geometry::MeshList MeshList;

        VersionId    version_id;
        std::fstream fs;
        std::string  directory;
        MeshList     meshes;
        bool         mesh_provided_as_interfaces;
        unsigned     nb_interfaces;
    };

    void GeomFile::load_meshes(Geometry& geometry) {
        
        fs.open(fname.c_str());

        if (!fs.is_open())
            throw OpenMEEG::OpenError(fname);

        //  Get the version of the geometry file format.

        unsigned major,minor; ///< version of the domain description
        fs >> io_utils::match("# Domain Description ") >> major >> io_utils::match(".") >> minor;

        if (fs.fail())
            throw OpenMEEG::WrongFileFormat(fname);

        version_id = version(major,minor);
        if (version_id==VERSION10)
            std::cerr << "(DEPRECATED) Please consider updating your geometry file to the new format 1.1 (see data/README.rst): "
                      << fname << std::endl;

        if (version_id==UNKNOWN_VERSION) {
             std::cerr << "Domain Description version not available !" << std::endl;
             throw OpenMEEG::WrongFileFormat(fname);
        }

        // Extract the absolute path of geometry file

        directory = absolute_path(fname);

        // Process meshes.

        bool has_meshfile;
        bool has_meshsection;

        if (version_id==VERSION11) {

            // Read the mesh section of the description file.
            // Try to load the meshfile (VTK::vtp file) or try to load the meshes

            fs >> io_utils::skip_comments("#") >> io_utils::match_optional("MeshFile",has_meshfile);

            if (has_meshfile) {
                GeometryIO* io = GeometryIO::create(filename());
                io->load(geometry);
            }

            fs >> io_utils::skip_comments("#") >> io_utils::match_optional("Meshes",has_meshsection);
            if (has_meshsection) {
                unsigned nb_meshes;
                fs >> nb_meshes;
                read_mesh_descriptions(nb_meshes,"Mesh");
            }
        }

        // Process interfaces.

        bool trash; // backward compatibility, catch "Mesh" optionally.
        fs >> io_utils::skip_comments('#')
            >> io_utils::match("Interfaces") >> nb_interfaces >> io_utils::match_optional("Mesh", trash);

        if (fs.fail())
            throw OpenMEEG::WrongFileFormat(fname);

        mesh_provided_as_interfaces = !has_meshfile && !has_meshsection;
        if (mesh_provided_as_interfaces)
            read_mesh_descriptions(nb_interfaces,"Interface");

        if (!has_meshfile)
            geometry.import(meshes);
    }

    void GeomFile::load_domains(Geometry& geometry) {

        // Create interfaces

        Interfaces interfaces;

        if (mesh_provided_as_interfaces) {
            for (const auto desc : meshes) {
                const std::string& name = desc.first;
                Interface interface(name);
                interface.oriented_meshes().push_back(OrientedMesh(geometry.mesh(name),OrientedMesh::Normal));
                interfaces.push_back(interface);
            }
        } else {
            for (unsigned i=0; i<nb_interfaces; ++i) {
                const std::string& interfacename = section_name(i,"Interface");
                interfaces.push_back(interfacename);
                for (auto& token : get_line_tokens()) {
                    const OrientedMesh::Orientation oriented = (extract_sign(token)) ? OrientedMesh::Normal : OrientedMesh::Opposite;
                    interfaces[i].oriented_meshes().push_back(OrientedMesh(geometry.mesh(token),oriented));
                }
            }
        }

        std::map<std::string,Interface&> interfaces_map;
        for (auto& interface : interfaces) {
            interfaces_map.insert({ interface.name(), interface });
            if (!interface.is_mesh_orientations_coherent()) { // check and correct global orientation
                std::cerr << "Interface \"" << interface.name() << "\" is not closed !" << std::endl
                          << "Please correct a mesh orientation when defining the interface in the geometry file." << std::endl;
                exit(1);
            }
        }

        //  Load domains.

        unsigned num_domains;
        fs >> io_utils::skip_comments('#') >> io_utils::match("Domains") >> num_domains;

        if (fs.fail())
            throw OpenMEEG::WrongFileFormat(fname);

        geometry.domains().resize(num_domains);
        for (auto& domain : geometry.domains()) {
            fs >> io_utils::skip_comments('#') >> io_utils::match("Domain");
            if (version_id==VERSION10) { // backward compatibility
                fs >> domain.name();
            } else {
                fs >> io_utils::token(domain.name(),':');
            }
            for (auto& token : get_line_tokens()) {
                if (token=="shared") {
                    std::cerr << "(DEPRECATED) Ignoring the useless shared keyword. Please consider updating the geometry file " << fname
                              << " to the new 1.1 format (see data/README.rst)." << std::endl;
                    break;                    
                }
                const SimpleDomain::Side side = (extract_sign(token)) ? SimpleDomain::Outside : SimpleDomain::Inside;
                try {
                    Interface& interface = interfaces_map.at(token);
                    domain.boundaries().push_back(SimpleDomain(interface,side));
                } catch(...) {
                    throw OpenMEEG::NonExistingDomain<std::string>(domain.name(),token);
                }
            }
        }

        if (fs.fail())
            throw OpenMEEG::WrongFileFormat(fname);
    }

    void GeomFile::save_geom(const Geometry& geometry) {
        fs.open(fname.c_str(),std::fstream::out);

        if (!fs.is_open())
            throw OpenMEEG::OpenError(fname);

        fs << "# Domain Description 1.1" << std::endl << std::endl;

        //  Create the list of interfaces.

        std::vector<const Mesh*> meshes;
        std::vector<const Interface*> interfaces;
        for (const auto& domain : geometry.domains())
            for (const auto& boundary : domain.boundaries()) {
                const Interface& interface = boundary.interface();
                if (std::find(interfaces.begin(),interfaces.end(),&interface)!=interfaces.end())
                    interfaces.push_back(&interface);
                    for (const auto& omesh : interface.oriented_meshes()) {
                        const Mesh& mesh = omesh.mesh();
                        if (std::find(meshes.begin(),meshes.end(),&mesh)!=meshes.end())
                            meshes.push_back(&mesh);
                    }
            }

        fs << "Meshes " << meshes.size() << std::endl << std::endl;
        
        for (const auto& mesh : meshes) {
            fs << "Mesh " << mesh->name() << ": \"" << "" << '"' << std::endl;
        }

        fs << "Interfaces " << interfaces.size() << std::endl << std::endl;

        for (const auto& interface : interfaces) {
            fs << "Interface " << interface->name() << ":";
            for (const auto& omesh : interface->oriented_meshes())
                fs << ' ' << ((omesh.orientation()==OrientedMesh::Normal) ? '+' : '-') << omesh.mesh().name();

            fs << std::endl;
        }
        fs << std::endl;

        fs << "Domains " << geometry.domains().size() << std::endl << std::endl;
        for (const auto& domain : geometry.domains()) {
            fs << "Domain " << domain.name() << ":";
            for (const auto& boundary : domain.boundaries())
                fs << ' ' << ((boundary.inside()) ? '-' : '+') << boundary.interface().name();
            fs << std::endl;
        }
        fs << std::endl;
    }
}
