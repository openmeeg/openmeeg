// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

#include <interface.h>
#include <OMExceptions.H>
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

        typedef enum { UNKNOWN_VERSION=-1, VERSION10, VERSION11 } VersionId;

        // Implementation of virtual methods.

        const char* name() const override { return "geom"; }

        void   load_meshes(Geometry& geometry) override;
        void   load_domains(Geometry& geometry) override;
        Matrix load_data() const override { return Matrix(); }

        virtual void save_geom(const Geometry& geometry) override;
        virtual void save_data(const Geometry&,const Matrix&) const override { }
        virtual void write() const override { }

        GeometryIO* clone(const std::string& filename) const override { return new GeomFile(filename); }

        // Helper functions.

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
            ifs >> io_utils::match(keyword) >> io_utils::token(res,':');
            return res;
        }

        std::string filename() {
            std::string filename;
            ifs >> io_utils::filename(filename,'"',false);
            return (std::filesystem::path(filename).is_relative()) ? (directory/filename).string() : filename;
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
            ifs >> io_utils::skip_comments("#") >> io_utils::match_optional(keyword+':',unnamed);
            const std::string& name = (unnamed || version_id==VERSION10) ? default_name(n+1) : given_name(keyword);
            return name;
        }

        Geometry::MeshList read_mesh_descriptions(const unsigned n,const std::string& keyword) {
            Geometry::MeshList meshes;
            for (unsigned i=0; i<n; ++i) {
                const std::string& name = section_name(i,keyword);
                const std::string& path = filename();
                meshes.push_back({ name, path });
            }
            return meshes;
        }

        //  Extract the tokens on the remainder of the line in the stream.

        std::vector<std::string> get_line_tokens() {
            std::string line;
            std::getline(ifs,line);
            std::istringstream iss(line);
            typedef std::istream_iterator<std::string> Iterator;
            std::vector<std::string> res;
            std::copy(Iterator(iss),Iterator(),std::back_inserter(res));
            return res;
        }

        GeomFile(const std::string& filename=""): base(filename,"geom") { }

        ~GeomFile() { ifs.close(); }

        static const GeomFile prototype;

        VersionId             version_id;
        std::ifstream         ifs;
        std::filesystem::path directory;
        bool                  mesh_provided_as_interfaces;
        unsigned              nb_interfaces;
    };

    void GeomFile::load_meshes(Geometry& geometry) {

        ifs.open(fname.c_str());

        if (!ifs.is_open())
            throw OpenMEEG::OpenError(fname);

        //  Get the version of the geometry file format.

        unsigned major,minor; ///< version of the domain description
        ifs >> io_utils::match("# Domain Description ") >> major >> io_utils::match(".") >> minor;

        if (ifs.fail())
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

        directory = std::filesystem::path(fname).parent_path();

        // Process meshes.

        bool has_meshfile = false;
        bool has_meshsection = false;

        if (version_id==VERSION11) {

            // Read the mesh section of the description file.
            // Try to load the meshfile (VTK::vtp file) or try to load the meshes

            ifs >> io_utils::skip_comments("#") >> io_utils::match_optional("MeshFile",has_meshfile);

            if (has_meshfile) {
                GeometryIO* io = GeometryIO::create(filename());
                io->load(geometry);
            }

            ifs >> io_utils::skip_comments("#") >> io_utils::match_optional("Meshes",has_meshsection);
            if (has_meshsection) {
                unsigned nb_meshes;
                ifs >> nb_meshes;
                const Geometry::MeshList& meshes = read_mesh_descriptions(nb_meshes,"Mesh");
                geometry.import(meshes);
            }
        }

        // Process interfaces.

        bool trash; // backward compatibility, catch "Mesh" optionally.
        ifs >> io_utils::skip_comments('#')
            >> io_utils::match("Interfaces") >> nb_interfaces >> io_utils::match_optional("Mesh", trash);

        if (ifs.fail())
            throw OpenMEEG::WrongFileFormat(fname);

        mesh_provided_as_interfaces = !has_meshfile && !has_meshsection;
        if (mesh_provided_as_interfaces) {
            const Geometry::MeshList& meshes = read_mesh_descriptions(nb_interfaces,"Interface");
            geometry.import(meshes);
        }
        // ideally we would close here, but load_domains comes immediately afterward, so leave it open
    }

    void GeomFile::load_domains(Geometry& geometry) {

        // Create interfaces

        Interfaces interfaces;

        if (mesh_provided_as_interfaces) {
            for (auto& mesh : geometry.meshes()) {
                const std::string& name = mesh.name();
                Interface interface(name);
                interface.oriented_meshes().push_back(OrientedMesh(mesh,OrientedMesh::Normal));
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
                throw OpenMEEG::WrongFileFormat(fname);
            }
        }

        //  Load domains.

        unsigned num_domains;
        ifs >> io_utils::skip_comments('#') >> io_utils::match("Domains") >> num_domains;

        if (ifs.fail())
            throw OpenMEEG::WrongFileFormat(fname);

        geometry.domains().resize(num_domains);
        for (auto& domain : geometry.domains()) {
            ifs >> io_utils::skip_comments('#') >> io_utils::match("Domain");
            if (version_id==VERSION10) { // backward compatibility
                ifs >> domain.name();
            } else {
                ifs >> io_utils::token(domain.name(),':');
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

        if (ifs.fail())
            throw OpenMEEG::WrongFileFormat(fname);

        ifs.close();
    }

    void GeomFile::save_geom(const Geometry& geometry) {
        std::ofstream ofs(fname.c_str());

        if (!ofs.is_open())
            throw OpenMEEG::OpenError(fname);

        ofs << "# Domain Description 1.1" << std::endl << std::endl;

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

        ofs << "Meshes " << meshes.size() << std::endl << std::endl;

        for (const auto& mesh : meshes)
            ofs << "Mesh " << mesh->name() << ": \"" << "" << '"' << std::endl;

        ofs << "Interfaces " << interfaces.size() << std::endl << std::endl;

        for (const auto& interface : interfaces) {
            ofs << "Interface " << interface->name() << ":";
            for (const auto& omesh : interface->oriented_meshes())
                ofs << ' ' << ((omesh.orientation()==OrientedMesh::Normal) ? '+' : '-') << omesh.mesh().name();

            ofs << std::endl;
        }
        ofs << std::endl;

        ofs << "Domains " << geometry.domains().size() << std::endl << std::endl;
        for (const auto& domain : geometry.domains()) {
            ofs << "Domain " << domain.name() << ":";
            for (const auto& boundary : domain.boundaries())
                ofs << ' ' << ((boundary.inside()) ? '-' : '+') << boundary.interface().name();
            ofs << std::endl;
        }
        ofs << std::endl;
    }
}
