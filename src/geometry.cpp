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

#include "om_utils.h"
#include "geometry.h"
#include "reader.h"
#include "PropertiesSpecialized.h"
#include <MeshDescription/Exceptions.H>
#include <IOUtils.H>


namespace OpenMEEG {

    const Domain Geometry::get_domain(const Vect3& p) const {

        for (Domains::const_iterator dit = this->domain_begin(); dit != this->domain_end(); dit++) {
            if (dit->contains_point(p)) {
                return *dit; // TODO 0 OR dit->name() OR *dit
            }
        }

    }

    void Geometry::read(const char* geomFileName, const char* condFileName) {

        destroy();
        has_cond() = false; // default parameter

        read_geom(geomFileName);

        // updates ? TODO
        geom_generate_indices();

        if(condFileName) {
            read_cond(condFileName);
            has_cond() = true;
        }
        has_cond() = true;
    }

    void Geometry::geom_generate_indices() { // this generates unique indices for vertices and triangles which will correspond to our unknowns.

        size_t index = 0;
        for ( Vertices::iterator pit = this->vertex_begin(); pit != this->vertex_end(); pit++, index++) {
            pit->index() = index;
        }
        // for ( Triangles::iterator tit = this->triangle_begin(); tit != this->triangle_end(); tit++, index++) { // TODO better triangle iterator ?
        for (iterator mit = this->begin(); mit != this->end(); mit++) {
            for (Mesh::iterator tit = mit->begin(); tit != mit->end(); tit++, index++) {
                tit->index() = index;
            }
        }
        this->size() = index;
    }

    void Geometry::read_cond(const char* condFileName) {

        typedef Utils::Properties::Named< std::string , Conductivity<double> > HeadProperties;
        HeadProperties properties(condFileName);

        // Store the internal conductivity of the external boundary of domain i
        // and store the external conductivity of the internal boundary of domain i
        for (Domains::iterator dit = this->domain_begin(); dit != domain_end(); dit++) {
            const Conductivity<double>& cond = properties.find(dit->name());
            dit->sigma() =  cond.sigma();
        }

        // std::cout << "\nChecking" << std::endl; 
        // for(int i = 0; i < nb_domains; i++) { // TODO print info on geom
            // std::cout << "\tMesh "  << " : internal conductivity = "  << sigin[i]
                // << " and external conductivity = " << sigout[i] << std::endl;
        // }
    }


    double Geometry::sigma(const std::string& name) const {
        for (std::vector<Domain>::const_iterator d = domain_begin(); d != domain_end(); d++) {
            if (name == d->name()) {
                return (d->sigma());
            }
        }
        return -1.; // TODO throw error unknownDomain
    }

    double Geometry::oriented(const Mesh& m1, const Mesh& m2) const {
        Domains doms = common_domains(m1, m2);
        double ans=0.;
        if (doms.size() == 2) { // TODO Maureen comment on the cylinder
            return 1.;
        }
        return (doms[0].meshOrient(m1) == doms[0].meshOrient(m2))?1.:-1.;
    }

    bool Geometry::selfCheck() const { // TODO: something else
        return true;
        // bool OK = true;
        // for(int i = 0; i < nb(); ++i)
        // {
            // const Mesh& m1 = getM(i);
            // if (!m1.has_correct_orientation()) {
                // warning(std::string("A mesh does not seem to be properly oriented"));
            // }
            // if(m1.has_self_intersection())
            // {
                // warning(std::string("Mesh is self intersecting !"));
                // m1.info();
                // OK = false;
                // std::cout << "Self intersection for mesh number " << i << std:: endl;
            // }
            // for(int j = i+1; j < nb(); ++j)
            // {
                // const Mesh& m2 = getM(j);
                // if(m1.intersection(m2))
                // {
                    // warning(std::string("2 meshes are intersecting !"));
                    // m1.info();
                    // m2.info();
                    // OK = false;
                // }
            // }
        // }
        // return OK;
    }

    /*
    // Mesh& Geometry::interior_mesh(const std::string& name) const {
        // for (std::vector<Domain>::const_iterator d = domain_begin(); d != domain_end(); d++) {
            // if ( name == d->name() ) {
                // return (this->begin()); // TODO Correct that now
            // }
        // }
        // return -1.; // TODO throw error unknownDomain
    // }


       bool Geometry::check(const Mesh& m) const {
       bool OK = true;
       if(m.has_self_intersection())
       {
       warning(std::string("Mesh is self intersecting !"));
       m.info();
       OK = false;
       }
       for(int i = 0; i < nb(); ++i)
       {
       const Mesh& m1 = getM(i);
       if(m1.intersection(m))
       {
       warning(std::string("Mesh is intersecting with one of the mesh in geom file !"));
       m1.info();
       OK = false;
       }
       }
       return OK;
       }

       int Geometry::getDomain(const Vect3& p) const {
       for (int i=0;i<nb();++i)
       if ((this->getM(i)).contains_vertex(p))
       return i;
       return nb();
       }
       */

    // For IO:s
    #ifdef USE_VTK
    void Geometry::get_data_from_vtk_reader(vtkPolyDataReader* reader, Mesh &m) {

        reader->Update();
        vtkPolyData *vtkMesh = reader->GetOutput();

        size_t npts;
        npts = vtkMesh->GetNumberOfPoints();

        vertices().reserve(npts + vertices().size()); // vertices
        normals().reserve(npts  + normals().size()); // normals on each point

        if (reader->GetNumberOfNormalsInFile()==0) {
            vtkPolyDataNormals *newNormals = vtkPolyDataNormals::New();
            newNormals->SetInput(vtkMesh);
            newNormals->Update();
            vtkMesh = newNormals->GetOutput();
        }

        vtkDataArray *normalsData = vtkMesh->GetPointData()->GetNormals();

        if (npts != normalsData->GetNumberOfTuples()) {
            std::cerr << "Error: number of vertices is not equal to number of normals in vtk file, correct or remove the normals." << std::endl;
            exit(1);
        }

        if (normalsData->GetNumberOfComponents() != 3) {
            std::cerr << "Error: wrong number of components of normals in vtk file, correct or remove the normals." << std::endl;
            exit(1);
        }

        for (int i = 0; i < npts; i++) {
            vertices().push_back(Vertex(vtkMesh->GetPoint(i)[0], vtkMesh->GetPoint(i)[1], vtkMesh->GetPoint(i)[2]));
            normals()[i](normalsData->GetTuple(i)[0], normalsData->GetTuple(i)[1], normalsData->GetTuple(i)[2]);
        }
        ntrgs = vtkMesh->GetNumberOfCells();
        m.reserve(ntrgs);
        vtkIdList *l;
        for (int i = 0; i<ntrgs; i++) {
            if (vtkMesh->GetCellType(i) == VTK_TRIANGLE) {
                l = vtkMesh->GetCell(i)->GetPointIds();
                m.push_back(Triangle(vertices()[l->GetId(0)+vertices().size()]-npts, vertices()[l->GetId(1)+vertices().size()-npts], vertices()[l->GetId(2)+vertices().size()-npts]))  ;
            } else {
                std::cerr << "This is not a triangulation" << std::endl;
                exit(1);
            }
        }
    }

    void Geometry::load_vtk(std::istream &is, Mesh &m) {
        destroy();

        // get length of file:
        is.seekg (0, ios::end);
        int length = is.tellg();
        is.seekg (0, ios::beg);

        // allocate memory:
        char * buffer = new char [length];

        // read data as a block:
        is.read (buffer, length);

        // held buffer by the array buf:
        vtkCharArray* buf = vtkCharArray::New();
        buf->SetArray(buffer, length, 1);

        vtkPolyDataReader* reader = vtkPolyDataReader::New();
        reader->SetInputArray(buf); // Specify 'buf' to be used when reading from a string
        reader->SetReadFromInputString(1);  // Enable reading from the InputArray 'buf' instead of the default, a file

        get_data_from_vtk_reader(reader);

        delete[] buffer;
        reader->Delete();

        m.update();
    }

    void Geometry::load_vtk(const char* filename, Mesh &m) {
        destroy();
        std::string s = filename;
        std::cout << "load_vtk : " << filename << std::endl;
        vtkPolyDataReader *reader = vtkPolyDataReader::New();
        reader->SetFileName(filename); // Specify file name of vtk data file to read
        if (!reader->IsFilePolyData()) {
            std::cerr << "This is not a valid vtk poly data file" << std::endl;
            reader->Delete();
            exit(1);
        }
        get_data_from_vtk_reader(reader);
        m.update();
    }

    void Geometry::load_vtp(std::istream &is) {
        destroy();
        // get length of file:
        is.seekg (0, ios::end);
        int length = is.tellg();
        is.seekg (0, ios::beg);

        // allocate memory:
        char * buffer = new char [length];

        // read data as a block:
        is.read (buffer, length);

        // held buffer by the array buf:
        vtkCharArray* buf = vtkCharArray::New();
        buf->SetArray(buffer, length, 1);

        vtkXMLPolyDataReader* reader = vtkXMLPolyDataReader::New();
        reader->SetInputArray(buf); // Specify 'buf' to be used when reading from a string
        reader->SetReadFromInputString(1);  // Enable reading from the InputArray 'buf' instead of the default, a file

        get_data_from_vtk_reader(reader);

        delete[] buffer;
        reader->Delete();

    }

    void Geometry::load_vtp(const char* filename) {
        destroy();
        std::string s = filename;
        std::cout << "load_vtp : " << filename << std::endl;
        vtkXMLPolyDataReader *reader = vtkXMLPolyDataReader::New();
        reader->SetFileName(filename); // Specify file name of vtk data file to read
        if (!reader->IsFilePolyData()) {
            std::cerr << "This is not a valid vtk poly data file" << std::endl;
            reader->Delete();
            exit(1);
        }

        get_data_from_vtk_reader(reader);
    }
    #endif

    #ifdef USE_GIFTI
    /*
    FIXME : implement GIFTI IOs
    */
    void from_gifti_image(gifti_image* gim) {
        std::cerr << "GIFTI reader : Not yet implemented" << std::endl;
        exit(1);
        return;
    }

    void Geometry::load_gifti(std::istream &is, Mesh &m) {
        std::cerr << "GIFTI reader : Not yet implemented" << std::endl;
        exit(1);
        destroy();
        // gifti_image* gim = gifti_read_image(filename);
        // from_gifti_image(gim)
        // update();
    }

    void Geometry::load_gifti(const char* filename, Mesh &m) {
        destroy();
        int read_data = 0;
        gifti_image* gim = gifti_read_image(filename, read_data);
        // from_gifti_image(gim);
        m.update();
        return;
    }
    #endif

    void Geometry::load_mesh_file(std::istream &is, Mesh &m) {
        destroy();

        unsigned char* uc = new unsigned char[5]; // File format
        is.read((char*)uc, sizeof(unsigned char)*5);
        delete[] uc;

        uc = new unsigned char[4]; // lbindian
        is.read((char*)uc, sizeof(unsigned char)*4);
        delete[] uc;

        unsigned int* ui = new unsigned int[1]; // arg_size
        is.read((char*)ui, sizeof(unsigned int));
        unsigned int arg_size = ui[0];
        delete[] ui;

        uc = new unsigned char[arg_size]; // Trash
        is.read((char*)uc, sizeof(unsigned char)*arg_size);
        delete[] uc;

        ui = new unsigned int[1]; // vertex_per_face
        is.read((char*)ui, sizeof(unsigned int));
        unsigned int vertex_per_face = ui[0];
        delete[] ui;

        ui = new unsigned int[1]; // mesh_time
        is.read((char*)ui, sizeof(unsigned int));
        unsigned int mesh_time = ui[0];
        delete[] ui;

        ui = new unsigned int[1]; // mesh_step
        is.read((char*)ui, sizeof(unsigned int));
        delete[] ui;

        ui = new unsigned int[1]; // vertex number
        is.read((char*)ui, sizeof(unsigned int));
        size_t npts = ui[0];
        delete[] ui;

        assert(vertex_per_face == 3); // Support only for triangulations
        assert(mesh_time == 1); // Support only 1 time frame

        float* pts_raw = new float[npts*3]; // Points
        is.read((char*)pts_raw, sizeof(float)*npts*3);

        ui = new unsigned int[1]; // arg_size
        is.read((char*)ui, sizeof(unsigned int));
        delete[] ui;

        float* normals_raw = new float[npts*3]; // Normals
        is.read((char*)normals_raw, sizeof(float)*npts*3);

        ui = new unsigned int[1]; // arg_size
        is.read((char*)ui, sizeof(unsigned int));
        delete[] ui;

        ui = new unsigned int[1]; // number of faces
        is.read((char*)ui, sizeof(unsigned int));
        size_t ntrgs;
        ntrgs = ui[0];
        delete[] ui;

        unsigned int* faces_raw = new unsigned int[ntrgs*3]; // Faces
        is.read((char*)faces_raw, sizeof(unsigned int)*ntrgs*3);

        vertices().reserve(npts + vertices().size()); // vertices
        normals().reserve(npts + normals().size()); // normals on each point
        for(int i = 0; i < m.nb_vertices(); ++i) {
            vertices().push_back(Vertex(pts_raw[i*3+0], pts_raw[i*3+1], pts_raw[i*3+2]));
            normals().push_back(Normal(normals_raw[i*3+0], normals_raw[i*3+1], normals_raw[i*3+2]));
        }
        m.reserve(ntrgs);
        for(int i = 0; i < ntrgs; ++i) {
            m.push_back(Triangle(vertices()[faces_raw[i*3+0]+vertices().size()-npts], vertices()[faces_raw[i*3+1]+vertices().size()-npts], vertices()[faces_raw[i*3+2]+vertices().size()-npts] ));
        }
        delete[] faces_raw;
        delete[] normals_raw;
        delete[] pts_raw;

        // update();
    }

    void Geometry::load_mesh_file(const char* filename, Mesh &m) {
        std::cout << "load_mesh_file : " << filename << std::endl;
        std::ifstream f(filename, std::ios::binary);
        if (!f.is_open()) {
            std::cerr << "Error opening MESH file: " << filename << std::endl;
            exit(1);
        }
        load_mesh_file(f, m);
        f.close();
    }

    void Geometry::load_tri_file(std::istream &f, Mesh &m) {

        destroy();

        f.seekg( 0, std::ios_base::beg );

        char ch;
        size_t npts, ntrgs;
        f >> ch;
        f >> npts;

        vertices().reserve(npts + vertices().size()); // vertices
        normals().reserve(npts + normals().size()); // normals on each point
        for (int i = 0; i < npts; i++) {
            Vertex p;
            Normal n;
            f >> p >> n;
            vertices().push_back(p);
            normals().push_back(n);
        }
        f >> ch >> ntrgs >> ntrgs >> ntrgs; // This number is repeated 3 times
        m.reserve(ntrgs);
        for (int i = 0; i < ntrgs; i++) {
            Triangle t;
            // f >> t;
            // TODO m[t];
        }

        m.update();
    }

    void Geometry::load_tri_file(const char* filename, Mesh &m) {

        std::string s = filename;
        std::cout << "load_tri : " << filename << std::endl;
        std::ifstream f(filename);
        if (!f.is_open()) {
            std::cerr << "Error opening TRI file: " << filename << std::endl;
            exit(1);
        }
        load_tri_file(f, m);
        f.close();
    }

    void Geometry::load_bnd_file(std::istream &f, Mesh &m) {

        destroy();

        std::string line;
        std::string st;

        f.seekg( 0, std::ios_base::beg );

        f >> io_utils::skip_comments('#') >> st;
        if (st == "Type=") {
            io_utils::skip_line(f);
            f >> io_utils::skip_comments('#') >> st;
        }

        assert(st == "NumberPositions=");
        size_t npts, ntrgs;
        f >> npts;

        f >> io_utils::skip_comments('#') >> st;
        if (st == "UnitPosition") {
            io_utils::skip_line(f); // skip : "UnitPosition mm"
        }

        f >> io_utils::skip_comments('#') >> st;
        assert(st == "Positions");

        vertices().reserve(npts + vertices().size()); // vertices

        for( int i = 0; i < npts; i++ ) {
            Vertex v;
            f >> io_utils::skip_comments('#') >> v;
            vertices().push_back(v);
        }

        f >> io_utils::skip_comments('#') >> st;
        assert(st == "NumberPolygons=");
        f >> io_utils::skip_comments('#') >> ntrgs;

        f >> io_utils::skip_comments('#') >> st;
        assert(st == "TypePolygons=");
        f >> io_utils::skip_comments('#') >> st;
        assert(st == "3");

        f >> io_utils::skip_comments('#') >> st;
        assert(st == "Polygons");

        m.reserve(ntrgs);
        for (int i=0; i<ntrgs; i++) {
             //TODO f >> io_utils::skip_comments('#') >> m[i];
        }

        m.update();

        normals().reserve(npts + normals().size()); // normals on each point
        // recompute_normals(); // Compute normals since bnd files don't have any !
    }

    void Geometry::load_bnd_file(const char* filename, Mesh &m) {
        std::string s = filename;
        std::cout << "load_bnd : " << filename << std::endl;
        std::ifstream f(filename);

        if (!f.is_open()) {
            std::cerr << "Error opening BND file: " << filename << std::endl;
            exit(1);
        }
        load_bnd_file(f, m);
        f.close();
    }

    void Geometry::load_off_file(std::istream &f, Mesh &m) {
        destroy();

        char tmp[128];
        int trash;
        f>>tmp;        // put the "OFF" string
        size_t npts, ntrgs;
        f>>npts;
        f>>ntrgs;
        f>>trash;

        vertices().reserve(npts + vertices().size()); // vertices
        for (int i = 0; i < npts; i++) {
            Vertex v;
            f >> v;
            vertices().push_back(v);
        }

        m.reserve(ntrgs);
        for (int i=0; i<ntrgs; i++) {
            f>>trash;        // put the "3" to trash
            // TODO f>>m[i];
        }

        m.update();

        normals().reserve(npts + normals().size()); // normals on each point
        // recompute_normals(); // Compute normals since off files don't have any !
    }

    void Geometry::load_off_file(const char* filename, Mesh &m) {

        std::string s = filename;
        std::cout << "load_off : " << filename << std::endl;
        std::ifstream f(filename);
        if (!f.is_open()) {
            std::cerr << "Error opening OFF file: " << filename << std::endl;
            exit(1);
        }
        load_off_file(f, m);
        f.close();
    }

    void Geometry::load_mesh(const char* filename, Mesh &m, const bool &verbose) {
        std::string extension = getNameExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
        if (extension == std::string("vtk")) {
            load_vtk_file(filename, m);
        } else if (extension == std::string("vtp")) {
            load_vtp_file(filename);
        } else if (extension == std::string("tri")) {
            load_tri_file(filename, m);
        } else if (extension == std::string("bnd")) {
            load_bnd_file(filename, m);
        } else if (extension == std::string("mesh")) {
            load_mesh_file(filename, m);
        } else if (extension == std::string("off")) {
            load_off_file(filename, m);
        } else if (extension == std::string("gii")) {
            load_gifti_file(filename, m);
        } else {
            std::cerr << "IO: load: Unknown mesh file format for " << filename << std::endl;
            exit(1);
        }

        if (verbose) {
            m.info();
        }
    }

    void Geometry::save_mesh(const char* filename, const Mesh& m) const {

        std::string extension = getNameExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
        if (extension==std::string("vtk")) {
            save_vtk_file(filename, m);
        } else if (extension==std::string("tri")) {
            save_tri_file(filename, m);
        } else if (extension==std::string("bnd")) {
            save_bnd_file(filename, m);
        } else if (extension==std::string("mesh")) {
            save_mesh_file(filename, m);
        } else if (extension==std::string("off")) {
            save_off_file(filename, m);
        } else {
            std::cerr << "Unknown file format for : " << filename << std::endl;
            exit(1);
        }
    }

    void Geometry::save_vtk_file(const char* filename, const Mesh& m) const {

        std::ofstream os(filename);
        os << "# vtk DataFile Version 2.0" << std::endl;
        os << "File " << filename << " generated by OpenMEEG" << std::endl;
        os << "ASCII" << std::endl;
        os << "DATASET POLYDATA" << std::endl;
        os << "POINTS " << m.nb_vertices() << " float" << std::endl;
        for(int i = 0; i < m.nb_vertices(); i++)
            os << m.vertices()[i] << std::endl;
        os << "POLYGONS " << m.nb_triangles() << " " << m.nb_triangles()*4 << std::endl;
        for(int i = 0; i < m.nb_triangles(); i++)
            os << 3 << " " << m[i] << std::endl;

        os << "CELL_DATA " << m.nb_triangles() << std::endl;
        os << "POINT_DATA " << m.nb_vertices() << std::endl;
        os << "NORMALS normals float" << std::endl;
        for(int i = 0; i < m.nb_vertices(); i++)
            os << normals()[i] << std::endl;

        os.close();
    }

    void Geometry::save_bnd_file(const char* filename, const Mesh& m) const {

        std::ofstream os(filename);
        os << "# Bnd mesh file generated by OpenMeeg" << std::endl;
        os << "Type= Unknown" << std::endl;
        os << "NumberPositions= " << m.nb_vertices() << std::endl;
        os << "UnitPosition\tmm" << std::endl;
        os << "Positions" << std::endl;
        for(int i = 0; i < m.nb_vertices(); i++)
            os << m.vertices()[i] << std::endl;

        os << "NumberPolygons= " << m.nb_triangles() << std::endl;
        os << "TypePolygons=\t3" << std::endl;
        os << "Polygons" << std::endl;
        for(int i=0; i<m.nb_triangles(); i++)
            os << m[i] << std::endl;

        os.close();
    }

    void Geometry::save_tri_file(const char* filename, const Mesh& m) const {

        std::ofstream os(filename);
        os << "- " << m.nb_vertices() << std::endl;
        for(int i=0; i<m.nb_vertices(); i++) {
            os << m.vertices()[i] << " ";
            os << normals()[i] << std::endl;
        }
        os << "- " << m.nb_triangles() << " " << m.nb_triangles() << " " << m.nb_triangles() << std::endl;
        for(int i=0; i<m.nb_triangles(); i++)
            os << m[i] << std::endl;

        os.close();
    }

    void Geometry::save_off_file(const char* filename, const Mesh& m) const {

        std::ofstream os(filename);
        os << "OFF" << std::endl;
        os << m.nb_vertices() << " " << m.nb_triangles() << " 0" << std::endl;
        for(int i=0; i<m.nb_vertices(); i++)
            os << m.vertices()[i] << std::endl;

        for(int i=0; i<m.nb_triangles(); i++)
            os << "3 " << m[i] << std::endl;

        os.close();
    }

    void Geometry::save_mesh_file(const char* filename, const Mesh& m) const {

        std::ofstream os(filename, std::ios::binary);

        unsigned char format[5] = {'b','i','n','a','r'}; // File format
        os.write((char*)format, sizeof(unsigned char)*5);

        unsigned char lbindian[4] = {'D','C','B','A'}; // lbindian
        os.write((char*)lbindian, sizeof(unsigned char)*4);

        unsigned int arg_size[1] = {4}; // arg_size
        os.write((char*)arg_size, sizeof(unsigned int));

        unsigned char VOID[4] = {'V','O','I','D'}; // Trash
        os.write((char*)VOID, sizeof(unsigned char)*4);

        unsigned int vertex_per_face[1] = {3}; // vertex_per_face
        os.write((char*)vertex_per_face, sizeof(unsigned int));

        unsigned int mesh_time[1] = {1}; // mesh_time
        os.write((char*)mesh_time, sizeof(unsigned int));

        unsigned int mesh_step[1] = {0}; // mesh_step
        os.write((char*)mesh_step, sizeof(unsigned int));

        unsigned int vertex_number[1] = {m.nb_vertices()}; // vertex number
        os.write((char*)vertex_number, sizeof(unsigned int));

        float* pts_raw = new float[m.nb_vertices()*3]; // Points
        float* normals_raw = new float[m.nb_vertices()*3]; // Normals
        unsigned int* faces_raw = new unsigned int[m.nb_triangles()*3]; // Faces

        for(int i = 0; i < m.nb_vertices(); ++i) {
            pts_raw[i*3+0] = (float)(m.vertices()[i]->x());
            pts_raw[i*3+1] = (float)(m.vertices()[i]->y());
            pts_raw[i*3+2] = (float)(m.vertices()[i]->z());
//            normals_raw[i*3+0] = (float)(m.normals()[i]->x());
//            normals_raw[i*3+1] = (float)(m.normals()[i]->y());
//            normals_raw[i*3+2] = (float)(m.normals()[i]->z());
        }
        for(int i = 0; i < m.nb_triangles(); ++i) {
 //           faces_raw[i*3+0] = m[i][0];
 //           faces_raw[i*3+1] = m[i][1];
 //           faces_raw[i*3+2] = m[i][2];
        }

        os.write((char*)pts_raw, sizeof(float)*m.nb_vertices()*3);          // vertices
        os.write((char*)vertex_number, sizeof(unsigned int));    // arg size : npts
        os.write((char*)normals_raw, sizeof(float)*m.nb_vertices()*3);      // normals
        unsigned char zero[1] = {0};
        os.write((char*)zero, sizeof(unsigned int));             // arg size : 0
        unsigned int faces_number[1] = {m.nb_triangles()};
        os.write((char*)faces_number, sizeof(unsigned int));     // ntrgs
        os.write((char*)faces_raw, sizeof(unsigned int)*m.nb_triangles()*3); // triangles

        delete[] faces_raw;
        delete[] normals_raw;
        delete[] pts_raw;
        os.close();
    }

    #ifdef USE_GIFTI
    gifti_image* to_gifti_image() {
        gifti_image* gim;
        return gim;
    }

    void Geometry::save_gifti_file(const char* filename, const Mesh& m) {
        std::cerr << "GIFTI writer : Not yet implemented" << std::endl;
        gifti_image* gim = to_gifti_image();
        int write_data = 1;
        if (gifti_write_image(gim, filename, write_data)) {
            std::cerr << "Error while writing GIFTI file" << std::endl;
            exit(1);
        }
        exit(1);
    }
    #endif
}
