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

namespace MeshDescription {

    #ifdef USE_VTK
    void get_data_from_vtk_reader(vtkPolyDataReader* reader) {

        reader->Update();
        vtkPolyData *vtkMesh = reader->GetOutput();

        npts = vtkMesh->GetNumberOfPoints();

        vertices.resize(npts + offset_vertices); // vertices
        normals.resize(npts+ offset_vertices); // normals on each point

        if (reader->GetNumberOfNormalsInFile()==0) {
            vtkPolyDataNormals *newNormals = vtkPolyDataNormals::New();
            newNormals->SetInput(vtkMesh);
            newNormals->Update();
            vtkMesh = newNormals->GetOutput();
        }

        vtkDataArray *normalsData = vtkMesh->GetPointData()->GetNormals();

        if (npts!=normalsData->GetNumberOfTuples()) {
            std::cerr << "Error: number of vertices is not equal to number of normals in vtk file, correct or remove the normals." << std::endl;
            exit(1);
        }

        if (normalsData->GetNumberOfComponents()!=3) {
            std::cerr << "Error: wrong number of components of normals in vtk file, correct or remove the normals." << std::endl;
            exit(1);
        }
            
        for (int i = 0; i<npts; i++) {
            vertices.push_back(vtkMesh->GetPoint(i)[0], vtkMesh->GetPoint(i)[1], vtkMesh->GetPoint(i)[2]);
            normals.push_back(normalsData->GetTuple(i)[0], normalsData->GetTuple(i)[1], normalsData->GetTuple(i)[2]);
        }
        ntrgs = vtkMesh->GetNumberOfCells();
        triangles.resize(ntrgs + offset_triangles);
        vtkIdList *l;
        for (int i = 0; i<ntrgs; i++) {
            if (vtkMesh->GetCellType(i) == VTK_TRIANGLE) {
                l = vtkMesh->GetCell(i)->GetPointIds();
                triangles.push_back(Triangle(vertices[l->GetId(0)+offset_vertices], vertices[l->GetId(1)+offset_vertices], vertices[l->GetId(2)+offset_vertices]))  ;
            } else {
                std::cerr << "This is not a triangulation" << std::endl;
                exit(1);
            }
        }
    }

    void load_vtk(std::istream &is) {
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

        update_triangles();
    }

    void load_vtk(const char* filename) {
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
        update_triangles();
    }

    void load_vtp(std::istream &is) {
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

        update_triangles();
    }

    void load_vtp(const char* filename) {
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
        update_triangles();
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

    void load_gifti(std::istream &is) {
        std::cerr << "GIFTI reader : Not yet implemented" << std::endl;
        exit(1);
        destroy();
        // gifti_image* gim = gifti_read_image(filename);
        // from_gifti_image(gim)
        // update_triangles();
    }

    void load_gifti(const char* filename) {
        destroy();
        int read_data = 0;
        gifti_image* gim = gifti_read_image(filename, read_data);
        // from_gifti_image(gim);
        update_triangles();
        return;
    }
    #endif

    void load_mesh(std::istream &is) {
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
        npts = ui[0];
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
        ntrgs = ui[0];
        delete[] ui;

        unsigned int* faces_raw = new unsigned int[ntrgs*3]; // Faces
        is.read((char*)faces_raw, sizeof(unsigned int)*ntrgs*3);

        vertices.resize(npts + offset_vertices); // vertices
        normals.resize(npts+ offset_vertices); // normals on each point
        for(int i = 0; i < npts; ++i) {
            vertices.push_back(pts_raw[i*3+0], pts_raw[i*3+1], pts_raw[i*3+2]);
            normals.push_back(normals_raw[i*3+0], normals_raw[i*3+1], normals_raw[i*3+2]);
        }
        triangles.resize(ntrgs + offset_triangles);
        for(int i = 0; i < ntrgs; ++i) {
            triangles[i+offset_triangles] = Triangle(vertices[faces_raw[i*3+0]+offset_vertices], vertices[face_raw[i*3+1]+offset_vertices], vertices[face_raw[i*3+2]+offset_vertices] );
        }
        delete[] faces_raw;
        delete[] normals_raw;
        delete[] pts_raw;

        // update_triangles();
    }

    void load_mesh(const char* filename) {
        std::cout << "load_mesh : " << filename << std::endl;
        std::ifstream f(filename,std::ios::binary);
        if (!f.is_open()) {
            std::cerr << "Error opening MESH file: " << filename << std::endl;
            exit(1);
        }
        load_mesh(f);
        f.close();
    }

    void load_tri(std::istream &f) {

        destroy();

        f.seekg( 0, std::ios_base::beg );

        char ch;
        f >> ch;
        f >> npts;

        vertices.resize(npts + offset_vertices); // vertices
        normals.resize(npts+ offset_vertices); // normals on each point
        for (int i = 0; i < npts; i++) {
            Point p;
            Normal n;
            f >> p >> n;
            vertices.push_back(p);
            normals.push_back(n);
        }
        f >> ch >> ntrgs >> ntrgs >> ntrgs; // This number is repeated 3 times
        triangles.resize(ntrgs + offset_triangles);
        for (int i = 0; i < ntrgs; i++) {
            Triangle t;
            f >> t;
            triangles.push_back(t);
        }

        update_triangles();
    }

    void load_tri(const char* filename) {

        std::string s = filename;
        std::cout << "load_tri : " << filename << std::endl;
        std::ifstream f(filename);
        if (!f.is_open()) {
            std::cerr << "Error opening TRI file: " << filename << std::endl;
            exit(1);
        }
        load_tri(f);
        f.close();
    }

    void load_bnd(std::istream &f) {

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
        f >> npts;

        f >> io_utils::skip_comments('#') >> st;
        if (st == "UnitPosition") {
            io_utils::skip_line(f); // skip : "UnitPosition mm"
        }

        f >> io_utils::skip_comments('#') >> st;
        assert(st == "Positions");

        vertices.resize(npts + offset_vertices); // vertices

        for( int i = 0; i < npts; i++ ) {
            f >> io_utils::skip_comments('#') >> vertices[i + offset_vertices];
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

        triangles.resize(ntrgs + offset_triangles);
        for (int i=0; i<ntrgs; i++) {
            f >> io_utils::skip_comments('#') >> triangles[i + offset_triangles];
        }

        update_triangles();

        normals.resize(npts + offset_vertices); // normals on each point
        recompute_normals(); // Compute normals since bnd files don't have any !
    }

    void load_bnd(const char* filename) {
        std::string s = filename;
        std::cout << "load_bnd : " << filename << std::endl;
        std::ifstream f(filename);

        if (!f.is_open()) {
            std::cerr << "Error opening BND file: " << filename << std::endl;
            exit(1);
        }
        load_bnd(f);
        f.close();
    }

    void load_off(std::istream &f) {
        destroy();

        char tmp[128];
        int trash;
        f>>tmp;        // put the "OFF" string
        f>>npts;
        f>>ntrgs;
        f>>trash;

        vertices.resize(npts + offset_vertices); // vertices
        for (int i=0; i<npts; i++) {
            f>>vertices[i + offset_vertices];
        }

        triangles.resize(ntrgs + offset_triangles);
        for (int i=0; i<ntrgs; i++) {
            f>>trash;        // put the "3" to trash
            f>>triangles[i + offset_triangles];
        }

        update_triangles();

        normals.resize(npts + offset_vertices); // normals on each point
        recompute_normals(); // Compute normals since off files don't have any !
    }

    void load_off(const char* filename) {

        std::string s = filename;
        std::cout << "load_off : " << filename << std::endl;
        std::ifstream f(filename);
        if (!f.is_open()) {
            std::cerr << "Error opening OFF file: " << filename << std::endl;
            exit(1);
        }
        load_off(f);
        f.close();
    }

    void save(const char* filename) const {

        std::string extension = getNameExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
        if (extension==std::string("vtk")) {
            save_vtk(filename);
        } else if (extension==std::string("tri")) {
            save_tri(filename);
        } else if (extension==std::string("bnd")) {
            save_bnd(filename);
        } else if (extension==std::string("mesh")) {
            save_mesh(filename);
        } else if (extension==std::string("off")) {
            save_off(filename);
        } else {
            std::cerr << "Unknown file format for : " << filename << std::endl;
            exit(1);
        }
    }

    void save_vtk(const char* filename) const {

        std::ofstream os(filename);
        os << "# vtk DataFile Version 2.0" << std::endl;
        os << "File " << filename << " generated by OpenMEEG" << std::endl;
        os << "ASCII" << std::endl;
        os << "DATASET POLYDATA" << std::endl;
        os << "POINTS " << npts << " float" << std::endl;
        for(int i=0; i<npts; i++)
            os << vertices[i] << std::endl;
        os << "POLYGONS " << ntrgs << " " << ntrgs*4 << std::endl;
        for(int i=0; i<ntrgs; i++)
            os << 3 << " " << triangles[i] << std::endl;

        os << "CELL_DATA " << ntrgs << std::endl;
        os << "POINT_DATA " << npts << std::endl;
        os << "NORMALS normals float" << std::endl;
        for(int i=0; i<npts; i++)
            os << normals[i] << std::endl;

        os.close();
    }

    void save_bnd(const char* filename) const {

        std::ofstream os(filename);
        os << "# Bnd mesh file generated by OpenMeeg" << std::endl;
        os << "Type= Unknown" << std::endl;
        os << "NumberPositions= " << npts << std::endl;
        os << "UnitPosition\tmm" << std::endl;
        os << "Positions" << std::endl;
        for(int i=0; i<npts; i++)
            os << vertices[i] << std::endl;

        os << "NumberPolygons= " << ntrgs << std::endl;
        os << "TypePolygons=\t3" << std::endl;
        os << "Polygons" << std::endl;
        for(int i=0; i<ntrgs; i++)
            os << triangles[i] << std::endl;

        os.close();
    }

    void save_tri(const char* filename) const {

        std::ofstream os(filename);
        os << "- " << npts << std::endl;
        for(int i=0; i<npts; i++) {
            os << vertices[i] << " ";
            os << normals[i] << std::endl;
        }
        os << "- " << ntrgs << " " << ntrgs << " " << ntrgs << std::endl;
        for(int i=0; i<ntrgs; i++)
            os << triangles[i] << std::endl;

        os.close();
    }

    void save_off(const char* filename) const {

        std::ofstream os(filename);
        os << "OFF" << std::endl;
        os << npts << " " << ntrgs << " 0" << std::endl;
        for(int i=0; i<npts; i++)
            os << vertices[i] << std::endl;

        for(int i=0; i<ntrgs; i++)
            os << "3 " << triangles[i] << std::endl;

        os.close();
    }

    void save_mesh(const char* filename) const {

        std::ofstream os(filename,std::ios::binary);

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

        unsigned int vertex_number[1] = {npts}; // vertex number
        os.write((char*)vertex_number, sizeof(unsigned int));

        float* pts_raw = new float[npts*3]; // Points
        float* normals_raw = new float[npts*3]; // Normals
        unsigned int* faces_raw = new unsigned int[ntrgs*3]; // Faces

        for(int i = 0; i < npts; ++i) {
            pts_raw[i*3+0] = (float)vertices[i].x();
            pts_raw[i*3+1] = (float)vertices[i].y();
            pts_raw[i*3+2] = (float)vertices[i].z();
            normals_raw[i*3+0] = (float)normals[i].x();
            normals_raw[i*3+1] = (float)normals[i].y();
            normals_raw[i*3+2] = (float)normals[i].z();
        }
        for(int i = 0; i < ntrgs; ++i) {
            faces_raw[i*3+0] = triangles[i][0];
            faces_raw[i*3+1] = triangles[i][1];
            faces_raw[i*3+2] = triangles[i][2];
        }

        os.write((char*)pts_raw, sizeof(float)*npts*3);          // vertices
        os.write((char*)vertex_number, sizeof(unsigned int));    // arg size : npts
        os.write((char*)normals_raw, sizeof(float)*npts*3);      // normals
        unsigned char zero[1] = {0};
        os.write((char*)zero, sizeof(unsigned int));             // arg size : 0
        unsigned int faces_number[1] = {ntrgs};
        os.write((char*)faces_number, sizeof(unsigned int));     // ntrgs
        os.write((char*)faces_raw, sizeof(unsigned int)*ntrgs*3); // triangles

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

    void save_gifti(const char* filename) {
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
