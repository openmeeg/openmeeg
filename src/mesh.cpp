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

#include "mesh.h"
#include "Triangle_triangle_intersection.h"

namespace OpenMEEG {

    // Mesh::Mesh(const Mesh& M) { *this = M; }

    std::istream& operator>>(std::istream &is, Mesh &m)
    {
        size_t a, b, c;
        is >> a >> b >> c;
        Triangle t(m.all_vertices()[a], m.all_vertices()[b], m.all_vertices()[c] );
        m.push_back(t);
        
        return is;
    }

    void Mesh::recompute_normals() {
        size_t i = 0;
        for (const_vertex_iterator vit = vertex_begin(); vit != vertex_end(); vit++, i++) {
            Normal normal(0);
            // for (SetTriangle::iterator tit = links[i].begin(); tit != links[i].end(); tit++) { TODO
            //     normal += tit->normal().normalize();
            // }
            normal.normalize();
            all_normals().push_back(normal);
        }
    }


    // Mesh& Mesh::operator= (const Mesh& M) {
        // if (this != &M)
            // copy(M);
        // return *this;
    // }

    // void Mesh::copy(const Mesh& M) {
        // destroy();
        // if (M.npts != 0) {
            // _all_vertices = M.all_vertices();
            // ntrgs   = M.ntrgs;
            // pts     = new Vect3[npts];
            // trgs    = new Triangle[ntrgs];
            // links   = new intSet[npts];
            // normals = new Vect3[npts];
            // for (int i=0; i < npts; i++) {
                // pts[i]     = M.pts[i];
                // links[i]   = M.links[i];
                // normals[i] = M.normals[i];
            // }
            // for (int i=0; i < ntrgs; i++)
                // trgs[i] = M.trgs[i];
        // }
    // }

    // void Mesh::destroy() {
        // if (npts != 0) {
            // delete [] pts;
            // delete [] trgs;
            // delete [] links;
            // delete [] normals;

            // npts = 0;
            // ntrgs = 0;
        // }
    // }

    void Mesh::update() {
        SetPVertex vset;
        for (Triangles::iterator tit = this->begin(); tit != this->end(); tit++) {
            tit->area()   = tit->normal().norm() / 2.0;
            vset.insert(&tit->s1().vertex());
            vset.insert(&tit->s2().vertex());
            vset.insert(&tit->s3().vertex());
        }
        // vertices().insert(vertex_begin(), vset.begin(), vset.end()); // copy the set of mesh vertices into a vector
        vset.clear();
        links().resize(vertices().size());
        size_t i = 0;
        for (const_vertex_iterator vit = vertex_begin(); vit != vertex_end(); vit++, i++) {
            for (const_iterator tit = this->begin(); tit != this->end(); tit++) {
                if (tit->contains(**vit)) {
                    links()[i].insert(*tit); 
                }
            }
        }
    }

    /**
     * Print informations about the mesh
    **/
    void Mesh::info() const {

        std::cout << "Mesh Info : "     << std::endl;
        std::cout << "\tName or ID : "  << name() << std::endl;
        std::cout << "\t# points : "    << nb_vertices() << std::endl;
        std::cout << "\t# triangles : " << nb_triangles() << std::endl;
        std::cout << "\tEuler characteristic : " << nb_vertices() - 3.*nb_triangles()/2. + nb_triangles() << std::endl;

        double min_area = std::numeric_limits<double>::max();
        double max_area = 0.;
        for (const_iterator tit = this->begin(); tit != this->end(); tit++) {
            min_area = std::min(tit->area(), min_area);
            max_area = std::max(tit->area(), max_area);
        }
        std::cout << "\tMin Area : " << min_area << std::endl;
        std::cout << "\tMax Area : " << max_area << std::endl;
    }

    // void Mesh::update_triangles() {
        // for (const_iterator tit = this->begin(); tit != this->end(); tit++) {
            // tit->normal() = (*tit)(0)->vertex().normal( (*tit)(1) , (*tit)(2) ); // TODO: not Normalized
            // tit->area()   = tit->normal().norm() / 2.0;
        // }
    // }

    bool Mesh::has_self_intersection() const {
        bool selfIntersects = false;
        for (const_iterator tit1 = this->begin(); tit1 != this->end(); tit1++) {
            for (const_iterator tit2 = tit1; tit2 != this->end(); tit2++) {
                if (!tit1->contains(tit2->s1().vertex()) && !tit1->contains(tit2->s2().vertex()) && !tit1->contains(tit1->s3().vertex())) {
                    if (triangle_intersection(*tit1, *tit2)) {
                        selfIntersects = true;
                        std::cout << "triangles " << tit1->index() << " and " << tit2->index() << " are intersecting" << std::endl;
                    }
                }
            }
        }
        return selfIntersects;
    }

    bool Mesh::intersection(const Mesh& m) const {
        bool intersects = false;
        for (const_iterator tit1 = this->begin(); tit1 != this->end(); tit1++) {
            for (const_iterator tit2 = m.begin(); tit2 != m.end(); tit2++) {
                intersects = intersects | triangle_intersection(*tit1, *tit2);
            }
        }
        return intersects;
    }

    bool Mesh::triangle_intersection(const Triangle& T1, const Triangle& T2 ) const {
        const Vect3& p1 = T1.s1().vertex();
        const Vect3& q1 = T1.s2().vertex();
        const Vect3& r1 = T1.s3().vertex();
        const Vect3& p2 = T2.s1().vertex();
        const Vect3& q2 = T2.s2().vertex();
        const Vect3& r2 = T2.s3().vertex();

        double pp1[3] = {p1.x(), p1.y(), p1.z()};
        double qq1[3] = {q1.x(), q1.y(), q1.z()};
        double rr1[3] = {r1.x(), r1.y(), r1.z()};
        double pp2[3] = {p2.x(), p2.y(), p2.z()};
        double qq2[3] = {q2.x(), q2.y(), q2.z()};
        double rr2[3] = {r2.x(), r2.y(), r2.z()};
        return tri_tri_overlap_test_3d(pp1, qq1, rr1, pp2, qq2, rr2);
    }

    // For IO:s
    #ifdef USE_VTK
    void Mesh::get_data_from_vtk_reader(vtkPolyDataReader* reader) {

        reader->Update();
        vtkPolyData *vtkMesh = reader->GetOutput();

        size_t npts;
        npts = vtkMesh->GetNumberOfPoints();

        all_vertices().reserve(npts + all_vertices().size()); // vertices
        all_normals().reserve(npts  + all_normals().size()); // normals on each point

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
            all_vertices().push_back(Vertex(vtkMesh->GetPoint(i)[0], vtkMesh->GetPoint(i)[1], vtkMesh->GetPoint(i)[2]));
            all_normals()[i](normalsData->GetTuple(i)[0], normalsData->GetTuple(i)[1], normalsData->GetTuple(i)[2]);
        }
        ntrgs = vtkMesh->GetNumberOfCells();
        m.reserve(ntrgs);
        vtkIdList *l;
        for (int i = 0; i<ntrgs; i++) {
            if (vtkMesh->GetCellType(i) == VTK_TRIANGLE) {
                l = vtkMesh->GetCell(i)->GetPointIds();
                m.push_back(Triangle(all_vertices()[l->GetId(0)+all_vertices().size()]-npts, all_vertices()[l->GetId(1)+all_vertices().size()-npts], all_vertices()[l->GetId(2)+all_vertices().size()-npts]))  ;
            } else {
                std::cerr << "This is not a triangulation" << std::endl;
                exit(1);
            }
        }
    }

    void Mesh::load_vtk_file(std::istream &is) {
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

    void Mesh::load_vtk(const char* filename) {
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

    void Mesh::load_vtp(std::istream &is) {
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

    void Mesh::load_vtp(const char* filename) {
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

    void Mesh::load_gifti(std::istream &is) {
        std::cerr << "GIFTI reader : Not yet implemented" << std::endl;
        exit(1);
        destroy();
        // gifti_image* gim = gifti_read_image(filename);
        // from_gifti_image(gim)
        // update();
    }

    void Mesh::load_gifti(const char* filename) {
        destroy();
        int read_data = 0;
        gifti_image* gim = gifti_read_image(filename, read_data);
        // from_gifti_image(gim);
        m.update();
        return;
    }
    #endif

    void Mesh::load_mesh_file(std::istream &is) {
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

        all_vertices().reserve(npts + all_vertices().size()); // vertices
        all_normals().reserve(npts + all_normals().size()); // normals on each point
        for(int i = 0; i < nb_vertices(); ++i) {
            add_vertex(Vertex(pts_raw[i*3+0], pts_raw[i*3+1], pts_raw[i*3+2]));
            add_normal(Normal(normals_raw[i*3+0], normals_raw[i*3+1], normals_raw[i*3+2]));
        }
        reserve(ntrgs);
        for(int i = 0; i < ntrgs; ++i) {
            push_back(Triangle(all_vertices()[faces_raw[i*3+0]+all_vertices().size()-npts], 
                                 all_vertices()[faces_raw[i*3+1]+all_vertices().size()-npts],
                                 all_vertices()[faces_raw[i*3+2]+all_vertices().size()-npts] ));
        }
        delete[] faces_raw;
        delete[] normals_raw;
        delete[] pts_raw;

        // update();
    }

    void Mesh::load_mesh_file(const char* filename) {
        std::cout << "load_mesh_file : " << filename << std::endl;
        std::ifstream f(filename, std::ios::binary);
        if (!f.is_open()) {
            std::cerr << "Error opening MESH file: " << filename << std::endl;
            exit(1);
        }
        load_mesh_file(f);
        f.close();
    }

    void Mesh::load_tri_file(std::istream &f) {

        destroy();

        f.seekg( 0, std::ios_base::beg );

        char ch;
        size_t npts, ntrgs;
        f >> ch;
        f >> npts;

        all_vertices().reserve(npts + all_vertices().size()); // vertices
        all_normals().reserve(npts + all_normals().size()); // normals on each point
        for (int i = 0; i < npts; i++) {
            Vertex v;
            Normal n;
            f >> v >> n;
            add_vertex(v);
            add_normal(n);
        }
        f >> ch >> ntrgs >> ntrgs >> ntrgs; // This number is repeated 3 times
        reserve(ntrgs);
        for (int i = 0; i < ntrgs; i++) {
            f >> *this;
        }

        update();
    }

    void Mesh::load_tri_file(const char* filename) {

        std::string s = filename;
        std::cout << "load_tri : " << filename << std::endl;
        std::ifstream f(filename);
        if (!f.is_open()) {
            std::cerr << "Error opening TRI file: " << filename << std::endl;
            exit(1);
        }
        load_tri_file(f);
        f.close();
    }

    void Mesh::load_bnd_file(std::istream &f) {

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

        all_vertices().reserve(npts + all_vertices().size()); // vertices

        for( int i = 0; i < npts; i++ ) {
            Vertex v;
            f >> io_utils::skip_comments('#') >> v;
            add_vertex(v);
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

        reserve(ntrgs);
        for ( size_t i = 0; i < ntrgs; i++) {
             f >> io_utils::skip_comments('#') >> *this;
        }

        update();

        all_normals().reserve(npts + all_normals().size()); // normals on each point
        recompute_normals(); // Compute normals since bnd files don't have any !
    }

    void Mesh::load_bnd_file(const char* filename) {
        std::string s = filename;
        std::cout << "load_bnd : " << filename << std::endl;
        std::ifstream f(filename);

        if (!f.is_open()) {
            std::cerr << "Error opening BND file: " << filename << std::endl;
            exit(1);
        }
        load_bnd_file(f);
        f.close();
    }

    void Mesh::load_off_file(std::istream &f) {
        destroy();

        char tmp[128];
        int trash;
        f>>tmp;        // put the "OFF" string
        size_t npts, ntrgs;
        f>>npts;
        f>>ntrgs;
        f>>trash;

        all_vertices().reserve(npts + all_vertices().size()); // vertices
        for (size_t i = 0; i < npts; i++) {
            Vertex v;
            f >> v;
            add_vertex(v);
        }

        reserve(ntrgs);
        for (size_t i = 0; i < ntrgs; i++) {
            f >> trash;        // put the "3" to trash
            f >> *this;
        }

        update();

        all_normals().reserve(npts + all_normals().size()); // normals on each point
        recompute_normals(); // Compute normals since off files don't have any !
    }

    void Mesh::load_off_file(const char* filename) {

        std::string s = filename;
        std::cout << "load_off : " << filename << std::endl;
        std::ifstream f(filename);
        if (!f.is_open()) {
            std::cerr << "Error opening OFF file: " << filename << std::endl;
            exit(1);
        }
        load_off_file(f);
        f.close();
    }

    void Mesh::load_mesh(const char* filename, const bool &verbose) {
        std::string extension = getNameExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
        if (extension == std::string("vtk")) {
            load_vtk_file(filename);
        } else if (extension == std::string("vtp")) {
            load_vtp_file(filename);
        } else if (extension == std::string("tri")) {
            load_tri_file(filename);
        } else if (extension == std::string("bnd")) {
            load_bnd_file(filename);
        } else if (extension == std::string("mesh")) {
            load_mesh_file(filename);
        } else if (extension == std::string("off")) {
            load_off_file(filename);
        } else if (extension == std::string("gii")) {
            load_gifti_file(filename);
        } else {
            std::cerr << "IO: load: Unknown mesh file format for " << filename << std::endl;
            exit(1);
        }

        if (verbose) {
            info();
        }
    }

    void Mesh::save_mesh(const char* filename) const {

        std::string extension = getNameExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
        if (extension==std::string("vtk")) {
            save_vtk_file(filename);
        } else if (extension==std::string("tri")) {
            save_tri_file(filename);
        } else if (extension==std::string("bnd")) {
            save_bnd_file(filename);
        } else if (extension==std::string("mesh")) {
            save_mesh_file(filename);
        } else if (extension==std::string("off")) {
            save_off_file(filename);
        } else {
            std::cerr << "Unknown file format for : " << filename << std::endl;
            exit(1);
        }
    }

    void Mesh::save_vtk_file(const char* filename) const {

        std::ofstream os(filename);
        os << "# vtk DataFile Version 2.0" << std::endl;
        os << "File " << filename << " generated by OpenMEEG" << std::endl;
        os << "ASCII" << std::endl;
        os << "DATASET POLYDATA" << std::endl;
        os << "POINTS " << nb_vertices() << " float" << std::endl;
        for(int i = 0; i < nb_vertices(); i++) {
            os << all_vertices()[i] << std::endl;
        }
        os << "POLYGONS " << nb_triangles() << " " << nb_triangles()*4 << std::endl;
        for(int i = 0; i < nb_triangles(); i++)
            os << 3 << " " << (*this)[i] << std::endl;

        os << "CELL_DATA " << nb_triangles() << std::endl;
        os << "POINT_DATA " << nb_vertices() << std::endl;
        os << "NORMALS normals float" << std::endl;
        for(int i = 0; i < nb_vertices(); i++)
            os << all_normals()[i] << std::endl;

        os.close();
    }

    void Mesh::save_bnd_file(const char* filename) const {

        std::ofstream os(filename);
        os << "# Bnd mesh file generated by OpenMeeg" << std::endl;
        os << "Type= Unknown" << std::endl;
        os << "NumberPositions= " << nb_vertices() << std::endl;
        os << "UnitPosition\tmm" << std::endl;
        os << "Positions" << std::endl;
        for(int i = 0; i < nb_vertices(); i++)
            os << all_vertices()[i] << std::endl;

        os << "NumberPolygons= " << nb_triangles() << std::endl;
        os << "TypePolygons=\t3" << std::endl;
        os << "Polygons" << std::endl;
        for(int i=0; i<nb_triangles(); i++)
            os << (*this)[i] << std::endl;

        os.close();
    }

    void Mesh::save_tri_file(const char* filename) const {

        std::ofstream os(filename);
        os << "- " << nb_vertices() << std::endl;
        for(int i=0; i<nb_vertices(); i++) {
            os << all_vertices()[i] << " ";
            os << all_normals()[i] << std::endl;
        }
        os << "- " << nb_triangles() << " " << nb_triangles() << " " << nb_triangles() << std::endl;
        for(int i=0; i<nb_triangles(); i++)
            os << (*this)[i] << std::endl;

        os.close();
    }

    void Mesh::save_off_file(const char* filename) const {

        std::ofstream os(filename);
        os << "OFF" << std::endl;
        os << nb_vertices() << " " << nb_triangles() << " 0" << std::endl;
        for(int i=0; i<nb_vertices(); i++)
            os << all_vertices()[i] << std::endl;

        for(int i=0; i<nb_triangles(); i++)
            os << "3 " << (*this)[i] << std::endl;

        os.close();
    }

    void Mesh::save_mesh_file(const char* filename) const {

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

        unsigned int vertex_number[1] = {nb_vertices()}; // vertex number
        os.write((char*)vertex_number, sizeof(unsigned int));

        float* pts_raw = new float[nb_vertices()*3]; // Points
        float* normals_raw = new float[nb_vertices()*3]; // Normals
        unsigned int* faces_raw = new unsigned int[nb_triangles()*3]; // Faces

        for(int i = 0; i < nb_vertices(); ++i) {
            pts_raw[i*3+0] = (float)(all_vertices()[i].x());
            pts_raw[i*3+1] = (float)(all_vertices()[i].y());
            pts_raw[i*3+2] = (float)(all_vertices()[i].z());
//            normals_raw[i*3+0] = (float)(all_normals()[i]->x());
//            normals_raw[i*3+1] = (float)(all_normals()[i]->y());
//            normals_raw[i*3+2] = (float)(all_normals()[i]->z());
        }
        for(int i = 0; i < nb_triangles(); ++i) {
 //           faces_raw[i*3+0] = *this[i][0];
 //           faces_raw[i*3+1] = *this[i][1];
 //           faces_raw[i*3+2] = *this[i][2];
        }

        os.write((char*)pts_raw, sizeof(float)*nb_vertices()*3);          // vertices
        os.write((char*)vertex_number, sizeof(unsigned int));    // arg size : npts
        os.write((char*)normals_raw, sizeof(float)*nb_vertices()*3);      // normals
        unsigned char zero[1] = {0};
        os.write((char*)zero, sizeof(unsigned int));             // arg size : 0
        unsigned int faces_number[1] = {nb_triangles()};
        os.write((char*)faces_number, sizeof(unsigned int));     // ntrgs
        os.write((char*)faces_raw, sizeof(unsigned int)*nb_triangles()*3); // triangles

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

    void Mesh::save_gifti_file(const char* filename) {
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
    /* TODO
    bool Mesh::has_correct_orientation() const {
        // First : check the local orientation (that all the triangles are all oriented in the same way)
        // define the triangle edges as (first point, second point)
        // if a triangle edge is ordered with (lower index, higher index) keep it in a edge_list
        // if not, exchange the two indices and put in a flipped_edge_list (lower index, higher index)
        // Transform the edge_list and flipped_edge_list unambigously into lists of numbers
        std::list<int> edge_list;
        std::list<int> flipped_edge_list;
        int radix = 10^(int(log2((double)npts)/log2(10.0))+1);
        for(int i=0;i<ntrgs;++i) {
            for(int j=1;j<=3;++j) {
                if (trgs[i].som(j)<trgs[i].next(j))
                    edge_list.push_back(trgs[i].som(j)*radix+trgs[i].next(j));
                else
                    flipped_edge_list.push_back(trgs[i].next(j)*radix+trgs[i].som(j));
            }
        }
        // Sort these two lists: they should be identical: if not, there is a local orientation problem.
        edge_list.sort();
        flipped_edge_list.sort();

        // check the global orientation: that the triangles are correctly oriented (outward-pointing normal)
        // compute the bounding box:
        double xmax = std::numeric_limits<int>::min();
        double ymax = std::numeric_limits<int>::min();
        double zmax = std::numeric_limits<int>::min();
        double xmin = std::numeric_limits<int>::max();
        double ymin = std::numeric_limits<int>::max();
        double zmin = std::numeric_limits<int>::max();
        for(int i=0; i<npts; ++i) {
            xmin = std::min(xmin, point(i).x());
            ymin = std::min(ymin, point(i).y());
            zmin = std::min(zmin, point(i).z());
            xmax = std::max(xmax, point(i).x());
            ymax = std::max(ymax, point(i).y());
            zmax = std::max(zmax, point(i).z());
        }
        Vect3 bbmin(xmin, ymin, zmin);
        Vect3 bbmax(xmax, ymax, zmax);
        Vect3 bbcenter = 0.5 * (bbmin + bbmax);

        // check the center of the bounding box is inside the mesh
        bool in_mesh = contains_point(bbcenter);
        if (!in_mesh)
            std::cerr << "Global orientation problem..." << std::endl << std::endl;

        if (flipped_edge_list != edge_list)
            std::cerr << "Local orientation problem..." << std::endl << std::endl;

        return in_mesh && (flipped_edge_list == edge_list);
    } */

} // end namespace OpenMeeg

