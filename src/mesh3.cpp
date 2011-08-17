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

#include <OpenMEEGConfigure.h>

#if WIN32
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_DEPRECATE 1
#endif
#include <math.h>
#include <assert.h>
#include "vect3.h"
#include "triangle.h"
#include "mesh3.h"
#include "om_utils.h"
#include "Triangle_triangle_intersection.h"
#include "IOUtils.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>

#ifdef USE_VTK
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataReader.h>
#include <vtkDataReader.h>
#include <vtkCellArray.h>
#include <vtkCharArray.h>
#include <vtkProperty.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkCleanPolyData.h>
#endif

#ifdef WIN32
inline double log2( double n )
{
    // log(n)/log(2) is log2.
    return log( n ) / log( 2.0 );
}
#endif

namespace OpenMEEG {

    Mesh::Mesh(): npts(0) {};

    Mesh::Mesh(const int a, const int b): npts(a), ntrgs(b), pts(new Vect3[npts]),
                            trgs(new Triangle[ntrgs]), links(new intSet[npts]),
                            normals(new Vect3[npts]) { }

    Mesh::Mesh(const Mesh& M): npts(0) { *this = M; }

    Mesh& Mesh::operator= (const Mesh& M) {
        if (this != &M)
            copy(M);
        return *this;
    }

    void Mesh::copy(const Mesh& M) {
        destroy();
        if (M.npts != 0) {
            npts    = M.npts;
            ntrgs   = M.ntrgs;
            pts     = new Vect3[npts];
            trgs    = new Triangle[ntrgs];
            links   = new intSet[npts];
            normals = new Vect3[npts];
            for (int i=0; i < npts; i++) {
                pts[i]     = M.pts[i];
                links[i]   = M.links[i];
                normals[i] = M.normals[i];
            }
            for (int i=0; i < ntrgs; i++)
                trgs[i] = M.trgs[i];
        }
    }

    void Mesh::destroy() {
        if (npts != 0) {
            delete [] pts;
            delete [] trgs;
            delete [] links;
            delete [] normals;

            npts = 0;
            ntrgs = 0;
        }
    }

    void Mesh::make_links() {
        for (int i=0; i<npts; i++)
            links[i].clear();

        for (int i=0; i<ntrgs; i++) {
            links[trgs[i].s1()].insert(i);
            links[trgs[i].s2()].insert(i);
            links[trgs[i].s3()].insert(i);
        }
    }

    void Mesh::get_file_format(const char* filename) {

        std::string extension = getNameExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
        if     (extension==std::string("vtk")) {
            streamFormat = Mesh::VTK;
        } else if(extension==std::string("tri")) {
            streamFormat = Mesh::TRI;
        } else if(extension==std::string("bnd")) {
            streamFormat = Mesh::BND;
        } else if(extension==std::string("mesh")) {
            streamFormat = Mesh::MESH;
        } else if(extension==std::string("off")) {
            streamFormat = Mesh::OFF;
        } else if(extension==std::string("gii")) {
            streamFormat = Mesh::GIFTI;
        } else {
            std::cerr << "Unknown mesh file format for " << filename << std::endl;
            exit(1);
        }
    }

    void Mesh::load(const char* filename, bool verbose) {

        std::string extension = getNameExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
        if (extension==std::string("vtk")) {
            load_vtk(filename);
        } else if(extension==std::string("tri")) {
            load_tri(filename);
        } else if(extension==std::string("bnd")) {
            load_bnd(filename);
        } else if(extension==std::string("mesh")) {
            load_mesh(filename);
        } else if(extension==std::string("off")) {
            load_off(filename);
        } else if(extension==std::string("gii")) {
            load_gifti(filename);
        } else {
            std::cerr << "Load : Unknown mesh file format for " << filename << std::endl;
            exit(1);
        }

        if(verbose)
            info();
    }

    #ifdef USE_VTK
    void Mesh::get_data_from_vtk_reader(vtkPolyDataReader* reader) {

        reader->Update();
        vtkPolyData *vtkMesh = reader->GetOutput();

        npts = vtkMesh->GetNumberOfPoints();

        pts = new Vect3[npts]; // points array
        normals = new Vect3[npts]; // normals on each point
        links = new intSet[npts];

        if(reader->GetNumberOfNormalsInFile()==0) {
            vtkPolyDataNormals *newNormals = vtkPolyDataNormals::New();
            newNormals->SetInput(vtkMesh);
            newNormals->Update();
            vtkMesh = newNormals->GetOutput();
        }

        vtkDataArray *normalsData = vtkMesh->GetPointData()->GetNormals();
        assert(npts == normalsData->GetNumberOfTuples());

        assert(3 == normalsData->GetNumberOfComponents());

        for (int i = 0; i<npts; i++) {
            pts[i].x() = vtkMesh->GetPoint(i)[0];
            pts[i].y() = vtkMesh->GetPoint(i)[1];
            pts[i].z() = vtkMesh->GetPoint(i)[2];

            normals[i](0) = normalsData->GetTuple(i)[0];
            normals[i](1) = normalsData->GetTuple(i)[1];
            normals[i](2) = normalsData->GetTuple(i)[2];
        }
        ntrgs = vtkMesh->GetNumberOfCells();
        trgs = new Triangle[ntrgs];
        vtkIdList *l;
        for (int i = 0; i<ntrgs; i++) {
            if (vtkMesh->GetCellType(i) == VTK_TRIANGLE) {
                l = vtkMesh->GetCell(i)->GetPointIds();
                trgs[i][0] = l->GetId(0);
                trgs[i][1] = l->GetId(1);
                trgs[i][2] = l->GetId(2);
            } else {
                std::cerr << "This is not a triangulation" << std::endl;
                exit(1);
            }
        }
    }

    void Mesh::load_vtk(std::istream &is) {
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

        Mesh::get_data_from_vtk_reader(reader);

        delete[] buffer;
        reader->Delete();

        make_links();
        update_triangles();
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

        Mesh::get_data_from_vtk_reader(reader);
        make_links();
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

    void Mesh::load_gifti(std::istream &is) {
        std::cerr << "GIFTI reader : Not yet implemented" << std::endl;
        exit(1);
        destroy();
        // gifti_image* gim = gifti_read_image(filename);
        // Mesh::from_gifti_image(gim)
        // make_links();
        // update_triangles();
    }

    void Mesh::load_gifti(const char* filename) {
        destroy();
        int read_data = 0;
        gifti_image* gim = gifti_read_image(filename, read_data);
        // Mesh::from_gifti_image(gim);
        make_links();
        update_triangles();
        return;
    }
    #endif

    void Mesh::load_mesh(std::istream &is) {
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

        pts = new Vect3[npts];
        normals = new Vect3[npts];
        links = new intSet[npts];
        for(int i = 0; i < npts; ++i) {
            pts[i].x() = pts_raw[i*3+0];
            pts[i].y() = pts_raw[i*3+1];
            pts[i].z() = pts_raw[i*3+2];
            normals[i].x() = normals_raw[i*3+0];
            normals[i].y() = normals_raw[i*3+1];
            normals[i].z() = normals_raw[i*3+2];
        }
        trgs = new Triangle[ntrgs];
        for(int i = 0; i < ntrgs; ++i) {
            trgs[i][0] = faces_raw[i*3+0];
            trgs[i][1] = faces_raw[i*3+1];
            trgs[i][2] = faces_raw[i*3+2];
        }
        delete[] faces_raw;
        delete[] normals_raw;
        delete[] pts_raw;

        make_links();
        update_triangles();
    }

    void Mesh::load_mesh(const char* filename) {
        std::string s = filename;
        std::cout << "load_mesh : " << filename << std::endl;
        std::ifstream f(filename);
        if(!f.is_open()) {
            std::cerr << "Error opening MESH file: " << filename << std::endl;
            exit(1);
        }
        Mesh::load_mesh(f);
        f.close();
    }

    void Mesh::load_tri(std::istream &f) {

        destroy();

        f.seekg( 0, std::ios_base::beg );

        char ch;
        f >> ch;
        f >> npts;

        pts = new Vect3[npts];
        normals = new Vect3[npts];
        links = new intSet[npts];
        for (int i=0; i<npts; i++) {
            f >> pts[i] >> normals[i];
        }
        f >> ch >> ntrgs >> ntrgs >> ntrgs; // This number is repeated 3 times
        trgs = new Triangle[ntrgs];
        for (int i=0; i<ntrgs; i++) {
            f >> trgs[i];
        }

        make_links();
        update_triangles();
    }

    void Mesh::load_tri(const char* filename) {

        std::string s = filename;
        std::cout << "load_tri : " << filename << std::endl;
        std::ifstream f(filename);
        if(!f.is_open()) {
            std::cerr << "Error opening TRI file: " << filename << std::endl;
            exit(1);
        }
        Mesh::load_tri(f);
        f.close();
    }

    void Mesh::load_bnd(std::istream &f) {

        destroy();

        std::string line;
        std::string st;

        f.seekg( 0, std::ios_base::beg );

        f >> io_utils::skip_comments('#') >> st;
        if(st == "Type=") {
            io_utils::skip_line(f);
            f >> io_utils::skip_comments('#') >> st;
        }

        assert(st == "NumberPositions=");
        f >> npts;

        f >> io_utils::skip_comments('#') >> st;
        if(st == "UnitPosition") {
            io_utils::skip_line(f); // skip : "UnitPosition mm"
        }

        f >> io_utils::skip_comments('#') >> st;
        assert(st == "Positions");

        pts = new Vect3[npts];
        links = new intSet[npts];

        for( int i = 0; i < npts; i += 1 ) {
            f >> io_utils::skip_comments('#') >> pts[i];
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

        trgs = new Triangle[ntrgs];
        for (int i=0; i<ntrgs; i++) {
            f >> io_utils::skip_comments('#') >> trgs[i];
        }

        make_links();
        update_triangles();

        normals = new Vect3[npts];
        recompute_normals(); // Compute normals since bnd files don't have any !
    }

    void Mesh::load_bnd(const char* filename) {
        std::string s = filename;
        std::cout << "load_bnd : " << filename << std::endl;
        std::ifstream f(filename);

        if(!f.is_open()) {
            std::cerr << "Error opening BND file: " << filename << std::endl;
            exit(1);
        }
        Mesh::load_bnd(f);
        f.close();
    }

    void Mesh::load_off(std::istream &f) {
        destroy();

        char tmp[128];
        int trash;
        f>>tmp;        // put the "OFF" string
        f>>npts;
        f>>ntrgs;
        f>>trash;

        pts = new Vect3[npts];
        links = new intSet[npts];
        for (int i=0; i<npts; i++) {
            f>>pts[i];
        }

        trgs = new Triangle[ntrgs];
        for (int i=0; i<ntrgs; i++) {
            f>>trash;        // put the "3" to trash
            f>>trgs[i];
        }

        make_links();
        update_triangles();

        normals = new Vect3[npts];
        recompute_normals(); // Compute normals since off files don't have any !
    }

    void Mesh::load_off(const char* filename) {

        std::string s = filename;
        std::cout << "load_off : " << filename << std::endl;
        std::ifstream f(filename);
        if(!f.is_open()) {
            std::cerr << "Error opening OFF file: " << filename << std::endl;
            exit(1);
        }
        Mesh::load_off(f);
        f.close();
    }

    void Mesh::save(const char* filename) const {

        std::string extension = getNameExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
        if     (extension==std::string("vtk")) {
            save_vtk(filename);
        } else if(extension==std::string("tri")) {
            save_tri(filename);
        } else if(extension==std::string("bnd")) {
            save_bnd(filename);
        } else if(extension==std::string("mesh")) {
            save_mesh(filename);
        } else if(extension==std::string("off")) {
            save_off(filename);
        } else {
            std::cerr << "Unknown file format for : " << filename << std::endl;
            exit(1);
        }
    }

    void Mesh::save_vtk(const char* filename) const {

        std::ofstream os(filename);
        os << "# vtk DataFile Version 2.0" << std::endl;
        os << "File " << filename << " generated by OpenMEEG" << std::endl;
        os << "ASCII" << std::endl;
        os << "DATASET POLYDATA" << std::endl;
        os << "POINTS " << npts << " float" << std::endl;
        for(int i=0; i<npts; i++)
            os << pts[i] << std::endl;
        os << "POLYGONS " << ntrgs << " " << ntrgs*4 << std::endl;
        for(int i=0; i<ntrgs; i++)
            os << 3 << " " << trgs[i] << std::endl;

        os << "CELL_DATA " << ntrgs << std::endl;
        os << "POINT_DATA " << npts << std::endl;
        os << "NORMALS normals float" << std::endl;
        for(int i=0; i<npts; i++)
            os << normals[i] << std::endl;

        os.close();
    }

    void Mesh::save_bnd(const char* filename) const {

        std::ofstream os(filename);
        os << "# Bnd mesh file generated by OpenMeeg" << std::endl;
        os << "Type= Unknown" << std::endl;
        os << "NumberPositions= " << npts << std::endl;
        os << "UnitPosition\tmm" << std::endl;
        os << "Positions" << std::endl;
        for(int i=0; i<npts; i++)
            os << pts[i] << std::endl;

        os << "NumberPolygons= " << ntrgs << std::endl;
        os << "TypePolygons=\t3" << std::endl;
        os << "Polygons" << std::endl;
        for(int i=0; i<ntrgs; i++)
            os << trgs[i] << std::endl;

        os.close();
    }

    void Mesh::save_tri(const char* filename) const {

        std::ofstream os(filename);
        os << "- " << npts << std::endl;
        for(int i=0; i<npts; i++) {
            os << pts[i] << " ";
            os << normals[i] << std::endl;
        }
        os << "- " << ntrgs << " " << ntrgs << " " << ntrgs << std::endl;
        for(int i=0; i<ntrgs; i++)
            os << trgs[i] << std::endl;

        os.close();
    }

    void Mesh::save_off(const char* filename) const {

        std::ofstream os(filename);
        os << "OFF" << std::endl;
        os << npts << " " << ntrgs << " 0" << std::endl;
        for(int i=0; i<npts; i++)
            os << pts[i] << std::endl;

        for(int i=0; i<ntrgs; i++)
            os << "3 " << trgs[i] << std::endl;

        os.close();
    }

    void Mesh::save_mesh(const char* filename) const {

        std::ofstream os(filename);

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
            pts_raw[i*3+0] = (float)pts[i].x();
            pts_raw[i*3+1] = (float)pts[i].y();
            pts_raw[i*3+2] = (float)pts[i].z();
            normals_raw[i*3+0] = (float)normals[i].x();
            normals_raw[i*3+1] = (float)normals[i].y();
            normals_raw[i*3+2] = (float)normals[i].z();
        }
        for(int i = 0; i < ntrgs; ++i) {
            faces_raw[i*3+0] = trgs[i][0];
            faces_raw[i*3+1] = trgs[i][1];
            faces_raw[i*3+2] = trgs[i][2];
        }

        os.write((char*)pts_raw, sizeof(float)*npts*3);          // points
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
    gifti_image* Mesh::to_gifti_image() {
        gifti_image* gim;
        return gim;
    }

    void Mesh::save_gifti(const char* filename) {
        std::cerr << "GIFTI writer : Not yet implemented" << std::endl;
        gifti_image* gim = Mesh::to_gifti_image();
        int write_data = 1;
        if (gifti_write_image(gim, filename, write_data)) {
            std::cerr << "Error while writing GIFTI file" << std::endl;
            exit(1);
        }
        exit(1);
    }
    #endif

    /**
     * Append another Mesh to on instance of Mesh
    **/
    void Mesh::append(const Mesh* m) {

        int old_npts = npts;
        int old_ntrgs = ntrgs;

        int new_npts = old_npts + m->nbPts();
        int new_ntrgs = old_ntrgs + m->nbTrgs();

        Vect3 *newPts = new Vect3[new_npts];
        Vect3 *newNormals = new Vect3[new_npts];
        Triangle *newTrgs = new Triangle[new_ntrgs];

        for(int i = 0; i < old_npts; ++i) {
            newPts[i] = point(i);
            newNormals[i] = normal(i);
        }

        for(int i = 0; i < m->nbPts(); ++i) {
            newPts[i + old_npts] = m->point(i);
            newNormals[i + old_npts] = m->normal(i);
        }

        for(int i = 0; i < old_ntrgs; ++i)
            newTrgs[i] = trgs[i];

        for(int i = 0; i < m->nbTrgs(); ++i) {
            const Triangle t = m->triangle(i);
            Triangle* new_t = new Triangle(t[0] + old_npts, t[1] + old_npts, t[2] + old_npts, t.normal());
            newTrgs[i + old_ntrgs] = *new_t;
        }

        destroy(); // Clean old Mesh
        npts = new_npts;
        ntrgs = new_ntrgs;
        pts = newPts;
        trgs = newTrgs;
        normals = newNormals;
        links = new intSet[npts];
        make_links(); // To keep a valid Mesh
        return;
    }

    /**
     * Get the neighboring triangle to triangle (a, b, c) containing edge (a, b)
    **/
    int Mesh::neighbor_triangle(int a, int b, int c) const {
        intSet possible_triangles = links[a];
        intSet::iterator it;
        for(it = possible_triangles.begin(); it != possible_triangles.end(); ++it) {
            Triangle t = triangle(*it);
            if (t.contains(b) && !t.contains(c)) {
                return *it;
            }
        }
        return -1; // Impossible to find neighboring triangle
    }

    /**
     * Print informations about the mesh
    **/
    void Mesh::info() const {

        std::cout << "Mesh Info : " << std::endl;
        std::cout << "\t# points : " << npts << std::endl;
        std::cout << "\t# triangles : " << ntrgs << std::endl;
        std::cout << "\tEuler characteristic : " << npts - 3*ntrgs/2 + ntrgs << std::endl;

        double min_area = trgs[0].area();
        double max_area = trgs[0].area();
        for(int i = 0; i < ntrgs; ++i) {
            min_area = std::min(trgs[i].area(),min_area);
            max_area = std::max(trgs[i].area(),max_area);
        }

        std::cout << "\tMin Area : " << min_area << std::endl;
        std::cout << "\tMax Area : " << max_area << std::endl;
    }

    /**
     * Smooth Mesh
    **/
    void Mesh::smooth(double smoothing_intensity, size_t niter) {

        std::vector< intSet > neighbors(npts);
        for(int i = 0; i < npts; ++i) {
            for(intSet::iterator it = links[i].begin(); it != links[i].end(); ++it) {
                Triangle nt = triangle(*it);
                for(size_t k = 0; k < 3; ++k) {
                    if (nt[k] != i) {
                        neighbors[i].insert(nt[k]);
                    }
                }
            }
        }
        Vect3* new_pts = new Vect3[npts];

        for(size_t n = 0; n < niter; ++n) {
            for(int p = 0; p < npts; ++p) {
                new_pts[p] = pts[p];
                for(intSet::iterator it = neighbors[p].begin(); it != neighbors[p].end(); ++it) {
                    new_pts[p] = new_pts[p] + (smoothing_intensity * (pts[*it] - pts[p])) / neighbors[p].size();
                }
            }
            for(int p = 0; p < npts; ++p) {
                pts[p] = new_pts[p];
            }
        }
        delete[] new_pts;
        update_triangles(); // Updating triangles (areas + normals)
        recompute_normals(); // Updating normals
    }

    /**
     * Surface Gradient
    **/
    SparseMatrix Mesh::gradient() const {

        SparseMatrix A(3*ntrgs, npts); // nb edges x points
        // loop on triangles
        for (int t=0; t<ntrgs; t++) {
            const Triangle& trg = triangle(t);
            Vect3 points[3] = {point(trg[0]), point(trg[1]), point(trg[2])};
            for(int j=0; j<3; j++) {
                Vect3 grads = P1Vector(points[0], points[1], points[2], j);
                for(int i=0; i<3; i++)
                    A(3*t+i,trg[j]) = grads(i);
            }
        }
        return A;
    }

    /**
     * Computes the total solid angle of a surface for a point p and tells whether p is inside the mesh or not.
     **/
    bool Mesh::contains_point(const Vect3& p) const {

        double solangle = 0.0;
        for (int itrg=0;itrg<ntrgs;itrg++)
            solangle += p.solangl(pts[trgs[itrg][0]],pts[trgs[itrg][1]],pts[trgs[itrg][2]]);

        if (std::abs(solangle) < 1e3*std::numeric_limits<double>::epsilon()) {
            return false;
        } else if (std::abs(solangle + 4*M_PI) < 1e3*std::numeric_limits<double>::epsilon()) {
            return true;
        } else if (std::abs(solangle - 4*M_PI) < 1e3*std::numeric_limits<double>::epsilon()) {
            std::cerr << "Mesh::contains_point(" << p << ") Error. Are you sure the mesh is properly oriented?\n";
            return false;
        } else {
            std::cerr << "Mesh::contains_point(" << p << ") Error. Are you sure the mesh is closed?\n"
                      << std::abs(solangle) << std::endl;
            exit(1);
        }
    }

    Vector Mesh::areas() const {

        Vector my_areas(nbTrgs());
        for(int i = 0; i < nbTrgs(); ++i)
            my_areas(i) = triangle(i).getArea();

        return my_areas;
    }

    void Mesh::update_triangles() {

        for(int i = 0; i < ntrgs; ++i) {
            trgs[i].normal() = pts[trgs[i][0]].normal( pts[trgs[i][1]] , pts[trgs[i][2]] );
            trgs[i].area() = trgs[i].normal().norm() / 2.0;
        }
    }

    void Mesh::recompute_normals() {
        for(int p = 0; p < nbPts(); ++p) {
            Vect3 my_normal(0);
            for(intSet::iterator it = links[p].begin(); it != links[p].end(); ++it)
                my_normal += trgs[*it].normal().normalize();

            my_normal.normalize();
            normals[p] = my_normal;
        }
    }

    bool Mesh::has_self_intersection() const {
        bool selfIntersects = false;
        for(int i = 0; i < ntrgs; ++i) {
            const Triangle& T1 = triangle(i);
            for(int j = i+1; j < ntrgs; ++j) {
                const Triangle& T2 = triangle(j);
                if (!T1.contains(T2.s1()) && !T1.contains(T2.s2()) && !T1.contains(T2.s3())) {
                    if (triangle_intersection(*this, i, *this, j)) {
                        selfIntersects = true;
                        std::cout << "triangles " << i << " and " << j << " are intersecting" <<std::endl;
                    }
                }
            }
        }
        return selfIntersects;
    }

    bool Mesh::intersection(const Mesh& m) const {
        bool intersects = false;
        for(int i = 0; i < ntrgs; ++i) {
            for(int j = 0; j < m.nbTrgs(); ++j) {
                intersects = intersects | triangle_intersection(*this, i, m, j);
            }
        }
        return intersects;
    }

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
    }

} // end namespace OpenMeeg
