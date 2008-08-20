/* FILE: $Id$ */

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
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

using namespace std;

Mesh::Mesh() {
    npts = 0;
}

Mesh::Mesh(int a, int b) {
    npts = a;
    ntrgs = b;
    pts = new Vect3[npts];
    trgs = new Triangle[ntrgs];
    links = new intSet[npts];
    normals = new Vect3[npts];
}

Mesh::Mesh(const Mesh& M) {
    *this = M;
}

Mesh& Mesh::operator=(const Mesh& M) {
    if (this!=&M) copy(M);
    return *this;
}

void Mesh::copy(const Mesh& M) {
    npts = M.npts;
    if (npts!=0) {
        ntrgs = M.ntrgs;
        pts = new Vect3[npts];
        trgs = new Triangle[ntrgs];
        links = new intSet[npts];
        normals = new Vect3[npts];
        for(int i=0;i<npts;i++){
            pts[i] = M.pts[i];
            links[i] = M.links[i];
        }
        for(int i=0;i<ntrgs;i++)
            trgs[i] = M.trgs[i];
    }
}

void Mesh::kill() {
    if (npts!=0) {
        delete []pts;
        delete []trgs;
        delete []links;
        delete []normals;

        npts = 0;
        ntrgs = 0;
    }
}

void Mesh::make_links() {
    int i;
    for (i=0;i<npts;i++) {
        links[i].clear();
    }

    for (i=0;i<ntrgs;i++) {
        links[trgs[i].s1()].insert(i);
        links[trgs[i].s2()].insert(i);
        links[trgs[i].s3()].insert(i);
    }
}

void Mesh::getFileFormat(const char* filename) {
    std::string extension = getNameExtension(filename);
    std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
    if     (extension==std::string("vtk"))  streamFormat = Mesh::VTK ;
    else if(extension==std::string("tri"))  streamFormat = Mesh::TRI;
    else if(extension==std::string("bnd"))  streamFormat = Mesh::BND;
    else if(extension==std::string("mesh")) streamFormat = Mesh::MESH;
    else if(extension==std::string("off"))  streamFormat = Mesh::OFF;
    else {
        cerr << "Unknown mesh file format for " << filename << endl;
        exit(1);
    }
}

void Mesh::load(const char* filename, bool checkClosedSurface, bool verbose) {
    std::string extension = getNameExtension(filename);
    std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
    if     (extension==std::string("vtk"))  load_vtk(filename,checkClosedSurface);
    else if(extension==std::string("tri"))  load_tri(filename,checkClosedSurface);
    else if(extension==std::string("bnd"))  load_bnd(filename,checkClosedSurface);
    else if(extension==std::string("mesh")) load_mesh(filename,checkClosedSurface);
    else if(extension==std::string("off"))  load_off(filename,checkClosedSurface);
    else {
        cerr << "Load : Unknown mesh file format for " << filename << endl;
        exit(1);
    }

    if(verbose)
    {
        info();
    }
}

#ifdef USE_VTK
void Mesh::getDataFromVTKReader(vtkPolyDataReader* reader) {   //private

    reader->Update();
    vtkPolyData *vtkMesh = reader->GetOutput();

    npts = vtkMesh->GetNumberOfPoints();

    pts = new Vect3[npts]; // points array
    normals = new Vect3[npts]; // normals on each point
    links = new intSet[npts];

    if(reader->GetNumberOfNormalsInFile()==0){
      vtkPolyDataNormals *newNormals = vtkPolyDataNormals::New();
      newNormals->SetInput(vtkMesh);
      newNormals->Update();
      vtkMesh = newNormals->GetOutput();
    }

    vtkDataArray *normalsData = vtkMesh->GetPointData()->GetNormals();
    assert(npts == normalsData->GetNumberOfTuples());

    assert(3 == normalsData->GetNumberOfComponents());

    for (int i = 0; i<npts; i++)
    {
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

void Mesh::load_vtk(std::istream &is, bool checkClosedSurface) {
    kill();

    // get length of file:
    is.seekg (0, ios::end);
    int length = is.tellg();
    is.seekg (0, ios::beg);

    // allocate memory:
    char * buffer = new char [length];

    // read data as a block:
    is.read (buffer,length);

    // held buffer by the array buf:
    vtkCharArray* buf = vtkCharArray::New();
    buf->SetArray(buffer,length,1);

    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetInputArray(buf); // Specify 'buf' to be used when reading from a string
    reader->SetReadFromInputString(1);  // Enable reading from the InputArray 'buf' instead of the default, a file

    Mesh::getDataFromVTKReader(reader);

    delete[] buffer;
    reader->Delete();

    make_links();
    update_triangles();
    updateTriangleOrientations(checkClosedSurface);
}

void Mesh::load_vtk(const char* filename, bool checkClosedSurface) {
    kill();

    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(filename); // Specify file name of vtk data file to read
    if (!reader->IsFilePolyData()) {
        std::cerr << "This is not a valid vtk poly data file" << std::endl;
        reader->Delete();
        exit(1);
    }

    Mesh::getDataFromVTKReader(reader);
}
#endif

void Mesh::load_mesh(std::istream &is, bool checkClosedSurface) {
    unsigned char* uc = new unsigned char[5]; // File format
    is.read((char*)uc,sizeof(unsigned char)*5);
    delete[] uc;

    uc = new unsigned char[4]; // lbindian
    is.read((char*)uc,sizeof(unsigned char)*4);
    delete[] uc;

    unsigned int* ui = new unsigned int[1]; // arg_size
    is.read((char*)ui,sizeof(unsigned int));
    unsigned int arg_size = ui[0];
    delete[] ui;

    uc = new unsigned char[arg_size]; // Trash
    is.read((char*)uc,sizeof(unsigned char)*arg_size);
    delete[] uc;

    ui = new unsigned int[1]; // vertex_per_face
    is.read((char*)ui,sizeof(unsigned int));
    unsigned int vertex_per_face = ui[0];
    delete[] ui;

    ui = new unsigned int[1]; // mesh_time
    is.read((char*)ui,sizeof(unsigned int));
    unsigned int mesh_time = ui[0];
    delete[] ui;

    ui = new unsigned int[1]; // mesh_step
    is.read((char*)ui,sizeof(unsigned int));
    delete[] ui;

    ui = new unsigned int[1]; // vertex number
    is.read((char*)ui,sizeof(unsigned int));
    npts = ui[0];
    delete[] ui;

    assert(vertex_per_face == 3); // Support only for triangulations
    assert(mesh_time == 1); // Support only 1 time frame

    float* pts_raw = new float[npts*3]; // Points
    is.read((char*)pts_raw,sizeof(float)*npts*3);

    ui = new unsigned int[1]; // arg_size
    is.read((char*)ui,sizeof(unsigned int));
    delete[] ui;

    float* normals_raw = new float[npts*3]; // Normals
    is.read((char*)normals_raw,sizeof(float)*npts*3);

    ui = new unsigned int[1]; // arg_size
    is.read((char*)ui,sizeof(unsigned int));
    delete[] ui;

    ui = new unsigned int[1]; // number of faces
    is.read((char*)ui,sizeof(unsigned int));
    ntrgs = ui[0];
    delete[] ui;

    unsigned int* faces_raw = new unsigned int[ntrgs*3]; // Faces
    is.read((char*)faces_raw,sizeof(unsigned int)*ntrgs*3);

    pts = new Vect3[npts];
    normals = new Vect3[npts];
    links = new intSet[npts];
    for(int i = 0; i < npts; ++i)
    {
        pts[i].x() = pts_raw[i*3+0];
        pts[i].y() = pts_raw[i*3+1];
        pts[i].z() = pts_raw[i*3+2];
        normals[i].x() = normals_raw[i*3+0];
        normals[i].y() = normals_raw[i*3+1];
        normals[i].z() = normals_raw[i*3+2];
    }
    trgs = new Triangle[ntrgs];
    for(int i = 0; i < ntrgs; ++i)
    {
        trgs[i][0] = faces_raw[i*3+0];
        trgs[i][1] = faces_raw[i*3+1];
        trgs[i][2] = faces_raw[i*3+2];
    }
    delete[] faces_raw;
    delete[] normals_raw;
    delete[] pts_raw;

    make_links();
    update_triangles();
    updateTriangleOrientations(checkClosedSurface);
}

void Mesh::load_mesh(const char* filename, bool checkClosedSurface) {
    string s = filename;
    cout << "load_mesh : " << filename << endl;
    std::ifstream f(filename);
    if(!f.is_open()) {
        cerr << "Error opening MESH file: " << filename << endl;
        exit(1);
    }
    Mesh::load_mesh(f,checkClosedSurface);
    f.close();
}

void Mesh::load_tri(std::istream &f, bool checkClosedSurface) {
    char ch; f>>ch;
    f>>npts;

    char myline[256];
    f.seekg( 0, std::ios_base::beg );
    f.getline(myline,256);
    f.getline(myline,256);
    int nread = 0;

    float r;
    nread = sscanf(myline,"%f %f %f %f %f %f",&r,&r,&r,&r,&r,&r);

    f.seekg( 0, std::ios_base::beg );
    f.getline(myline,256);

    pts = new Vect3[npts];
    normals = new Vect3[npts];
    links = new intSet[npts];
    for (int i=0;i<npts;i++){
        f>>pts[i];
        pts[i] = pts[i];
        f >> normals[i];
    }
    f>>ch;
    f>>ntrgs; f>>ntrgs; f>>ntrgs; // This number is repeated 3 times
    trgs = new Triangle[ntrgs];
    for (int i=0;i<ntrgs;i++) {
        f>>trgs[i];
    }

    make_links();
    update_triangles();
    updateTriangleOrientations(checkClosedSurface);
}

void Mesh::load_tri(const char* filename, bool checkClosedSurface) {
    string s = filename;
    cout << "load_tri : " << filename << endl;
    std::ifstream f(filename);
    if(!f.is_open()) {
        cerr << "Error opening TRI file: " << filename << endl;
        exit(1);
    }
    Mesh::load_tri(f,checkClosedSurface);
    f.close();
}



void Mesh::load_bnd(std::istream &f, bool checkClosedSurface) {
    std::string line;
    string st;

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
    for (int i=0;i<ntrgs;i++) {
        f >> io_utils::skip_comments('#') >> trgs[i];
    }

    make_links();
    update_triangles();
    updateTriangleOrientations(checkClosedSurface);

    normals = new Vect3[npts];
    recompute_normals(); // Compute normals since bnd files don't have any !
}

void Mesh::load_bnd(const char* filename, bool checkClosedSurface) {
    string s = filename;
    cout << "load_bnd : " << filename << endl;
    std::ifstream f(filename);

    if(!f.is_open()) {
        cerr << "Error opening BND file: " << filename << endl;
        exit(1);
    }
    Mesh::load_bnd(f);
    f.close();
}

void Mesh::load_off(std::istream &f, bool checkClosedSurface) {

    kill();
    char tmp[128]; int trash;
    f>>tmp;        // put the "OFF" string
    f>>npts; f>>ntrgs; f>>trash;

    pts = new Vect3[npts];
    links = new intSet[npts];
    for (int i=0;i<npts;i++) {
        f>>pts[i];
    }

    trgs = new Triangle[ntrgs];
    for (int i=0;i<ntrgs;i++) {
        f>>trash;        // put the "3" to trash
        f>>trgs[i];
    }

    make_links();
    update_triangles();
    updateTriangleOrientations(checkClosedSurface);

    normals = new Vect3[npts];
    recompute_normals(); // Compute normals since off files don't have any !
}

void Mesh::load_off(const char* filename, bool checkClosedSurface) {
    string s = filename;
    cout << "load_off : " << filename << endl;
    std::ifstream f(filename);
    if(!f.is_open()) {
        cerr << "Error opening OFF file: " << filename << endl;
        exit(1);
    }
    Mesh::load_off(f);
    f.close();
}

void Mesh::save(const char* filename) {
    std::string extension = getNameExtension(filename);
    std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
    if     (extension==std::string("vtk"))  save_vtk(filename);
    else if(extension==std::string("tri"))  save_tri(filename);
    else if(extension==std::string("bnd"))  save_bnd(filename);
    else if(extension==std::string("mesh")) save_mesh(filename);
    else if(extension==std::string("off"))  save_off(filename);
    else {
        cerr << "Unknown file format for : " << filename << endl;
        exit(1);
    }
}

void Mesh::save_vtk(const char* filename) {
    std::ofstream os(filename);
    os<<"# vtk DataFile Version 2.0"<<std::endl;
    os<<"File "<<filename<<" generated by OpenMEEG"<<std::endl;
    os<<"ASCII"<<std::endl;
    os<<"DATASET POLYDATA"<<std::endl;
    os<<"POINTS "<<npts<<" float"<<std::endl;
    for(int i=0;i<npts;i++) os<<pts[i]<<std::endl;
    os<<"POLYGONS "<<ntrgs<<" "<<ntrgs*4<<std::endl;
    for(int i=0;i<ntrgs;i++) os<<3<<" "<<trgs[i]<<std::endl;

    os<<"CELL_DATA "<<ntrgs<<std::endl;
    os<<"POINT_DATA "<<npts<<std::endl;
    os<<"NORMALS normals float"<<std::endl;
    for(int i=0;i<npts;i++)
        os<<normals[i]<<std::endl;
    os.close();
}

void Mesh::save_bnd(const char* filename) {
    std::ofstream os(filename);
    os << "# Bnd mesh file generated by OpenMeeg" << std::endl;
    os << "Type= Unknown" << std::endl;
    os << "NumberPositions= " << npts << std::endl;
    os << "UnitPosition\tmm" << std::endl;
    os << "Positions" << std::endl;
    for(int i=0;i<npts;i++) os<<pts[i]<<std::endl;
    os << "NumberPolygons= " << ntrgs <<std::endl;
    os << "TypePolygons=\t3" << std::endl;
    os << "Polygons" << std::endl;
    for(int i=0;i<ntrgs;i++) os << trgs[i] << std::endl;
    os.close();
}

void Mesh::save_tri(const char* filename) {
    std::ofstream os(filename);
    os<<"- "<<this->npts<<std::endl;
    for(int i=0;i<npts;i++) {
        os<<this->pts[i]<<" ";
        os<<this->normals[i]<<std::endl;
    }
    os<<"- "<<ntrgs<<" "<<ntrgs<<" "<<ntrgs<<std::endl;
    for(int i=0;i<ntrgs;i++) os<<this->trgs[i]<<std::endl;
    os.close();
}

void Mesh::save_off(const char* filename) {
    std::ofstream os(filename);
    os<<"OFF"<<std::endl;
    os<<this->npts<<" "<<this->ntrgs<<" 0"<<std::endl;
    for(int i=0;i<npts;i++) {
        os<<this->pts[i]<<std::endl;
    }
    for(int i=0;i<ntrgs;i++) os<<"3 "<<this->trgs[i]<<std::endl;
    os.close();
}

void Mesh::save_mesh(const char* filename) {
    std::ofstream os(filename);

    unsigned char format[5] = {'b','i','n','a','r'}; // File format
    os.write((char*)format,sizeof(unsigned char)*5);

    unsigned char lbindian[4] = {'D','C','B','A'}; // lbindian
    os.write((char*)lbindian,sizeof(unsigned char)*4);

    unsigned int arg_size[1] = {4}; // arg_size
    os.write((char*)arg_size,sizeof(unsigned int));

    unsigned char VOID[4] = {'V','O','I','D'}; // Trash
    os.write((char*)VOID,sizeof(unsigned char)*4);

    unsigned int vertex_per_face[1] = {3}; // vertex_per_face
    os.write((char*)vertex_per_face,sizeof(unsigned int));

    unsigned int mesh_time[1] = {1}; // mesh_time
    os.write((char*)mesh_time,sizeof(unsigned int));

    unsigned int mesh_step[1] = {0}; // mesh_step
    os.write((char*)mesh_step,sizeof(unsigned int));

    unsigned int vertex_number[1] = {npts}; // vertex number
    os.write((char*)vertex_number,sizeof(unsigned int));

    float* pts_raw = new float[npts*3]; // Points
    float* normals_raw = new float[npts*3]; // Normals
    unsigned int* faces_raw = new unsigned int[ntrgs*3]; // Faces

    for(int i = 0; i < npts; ++i)
    {
        pts_raw[i*3+0] = (float)pts[i].x();
        pts_raw[i*3+1] = (float)pts[i].y();
        pts_raw[i*3+2] = (float)pts[i].z();
        normals_raw[i*3+0] = (float)normals[i].x();
        normals_raw[i*3+1] = (float)normals[i].y();
        normals_raw[i*3+2] = (float)normals[i].z();
    }
    for(int i = 0; i < ntrgs; ++i)
    {
        faces_raw[i*3+0] = trgs[i][0];
        faces_raw[i*3+1] = trgs[i][1];
        faces_raw[i*3+2] = trgs[i][2];
    }

    os.write((char*)pts_raw,sizeof(float)*npts*3);          // points
    os.write((char*)vertex_number,sizeof(unsigned int));    // arg size : npts
    os.write((char*)normals_raw,sizeof(float)*npts*3);      // normals
    unsigned char zero[1] = {0};
    os.write((char*)zero,sizeof(unsigned int));             // arg size : 0
    unsigned int faces_number[1] = {ntrgs};
    os.write((char*)faces_number,sizeof(unsigned int));     // ntrgs
    os.write((char*)faces_raw,sizeof(unsigned int)*ntrgs*3); // triangles

    delete[] faces_raw;
    delete[] normals_raw;
    delete[] pts_raw;
    os.close();
}

/**
 * Append another Mesh to on instance of Mesh
**/
void Mesh::append(const Mesh* m) {
    int old_npts = this->npts;
    int old_ntrgs = this->ntrgs;

    int new_npts = old_npts + m->nbPts();
    int new_ntrgs = old_ntrgs + m->nbTrgs();

    Vect3 *newPts = new Vect3[new_npts];
    Vect3 *newNormals = new Vect3[new_npts];
    Triangle *newTrgs = new Triangle[new_ntrgs];

    for(int i = 0; i < old_npts; ++i) {
        newPts[i] = this->getPt(i);
        newNormals[i] = this->normal(i);
    }

    for(int i = 0; i < m->nbPts(); ++i) {
        newPts[i + old_npts] = m->getPt(i);
        newNormals[i + old_npts] = m->normal(i);
    }

    for(int i = 0; i < old_ntrgs; ++i)
        newTrgs[i] = this->trgs[i];

    for(int i = 0; i < m->nbTrgs(); ++i) {
        const Triangle t = m->getTrg(i);
        Triangle* new_t = new Triangle(t[0] + old_npts, t[1] + old_npts, t[2] + old_npts, t.normal());
        newTrgs[i + old_ntrgs] = *new_t;
    }

    kill(); // Clean old Mesh
    npts = new_npts;
    ntrgs = new_ntrgs;
    pts = newPts;
    trgs = newTrgs;
    normals = newNormals;
    links = new intSet[npts];
    make_links(); // To keep a valid Mesh
    //updateTriangleOrientations(true);
    return;
}

/**
 * Update the orientations of the triangles in the mesh to have all
 * normals pointing in the same direction
**/
void Mesh::updateTriangleOrientations(bool checkClosedSurface) {
    return;
    // std::cout << "Updating Triangle Orientations" << std::endl;
    std::vector< bool > seen(ntrgs,false); // Flag to say if a Triangle has been seen or not
    std::list< int > triangles;
    triangles.push_front(0); // Add First Triangle to the Heap
    seen[0] = true;

    ncomponents = 1;
    ninversions = 0;

    for(int i = 0; i < ntrgs; ++i)
    {
        if(triangles.empty()) // Try to find an unseen triangle in an other connexe component
        {
            int j;
            for(j = 0; j < ntrgs; ++j)
            {
                if(!seen[j])
                {
                    triangles.push_front(j);
                    seen[j] = true;
                    ncomponents++;
                    break;
                }
            }
            if(j == ntrgs)
            {
                std::cerr << "Problem while updating triangles orientations !"<< std::endl;
                exit(1);
            }
        }

        int next = triangles.front();
        triangles.pop_front();
        Triangle t = getTrg(next);

        for(int offset = 0; offset < 3; ++offset)
        {
            int a = t.next(offset);
            int b = t.next(offset+1);
            int c = t.next(offset+2);
            int n = getNeighTrg(a,b,c);
            if(n >= 0)
            {
                Triangle nt = getTrg(n);

                if(nt.getArea() < 0.000001) // TODO : check what should be a good value for the threshold
                {
                    std::cerr << "Error : Mesh contains a flat triangle !" << std::endl;
                    exit(1);
                }

                int aId = trgs[n].contains(a);
                int bId = trgs[n].contains(b);
                int dId;
                for(int k = 1; k <= 3; ++k)
                {
                    if ( (aId != k) && (bId != k)) dId = k;
                }
                int d = trgs[n].som(dId);
                trgs[n].s1() = b;
                trgs[n].s2() = a;
                trgs[n].s3() = d;

                // Check inversion :
                Vect3 old_normal = (getPt(nt.s2())-getPt(nt.s1()))^(getPt(nt.s3())-getPt(nt.s1()));
                Vect3 new_normal = (getPt(trgs[n].s2())-getPt(trgs[n].s1()))^(getPt(trgs[n].s3())-getPt(trgs[n].s1()));
                if((old_normal * new_normal) <= 0)
                {
                    ninversions++;
                    std::cout << "Warning : Fixing triangle " << n << std::endl;
                    std::cout << "\tOld triangle : " << nt << std::endl;
                    std::cout << "\tNew triangle : " << trgs[n] << std::endl;
                }

                if(!seen[n]) {
                    triangles.push_front(n);
                    seen[n] = true;
                }
            } else {
                if(checkClosedSurface)
                {
                    std::cerr << "!!!!!!!!!!!! Warning : Impossible to find neighboring triangle !!!!!!!!!!!!" << std::endl;
                    std::cerr << "!!!!!!!!!!!! Your mesh does not seem to be a closed surface    !!!!!!!!!!!!" << std::endl;
                    exit(1);
                }
            }
        }
    }
    for(int i = 0; i < ntrgs; ++i) assert(seen[i]); // All triangles have been explored
    assert(triangles.empty());
    return;
}

/**
 * Get the neighboring triangle to triangle (a,b,c) containing edge (a,b)
**/
int Mesh::getNeighTrg(int a, int b, int c) const {
    intSet possible_triangles = links[a];
    intSet::iterator it;
    for(it = possible_triangles.begin(); it != possible_triangles.end(); ++it) {
        Triangle t = getTrg(*it);
        if (t.contains(b) && !t.contains(c)) return *it;
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
    // std::cout << "\tNb of connexe components : " << ncomponents << std::endl;
    // std::cout << "\tNb of Triangles inversions : " << ninversions << std::endl;

    double min_area = trgs[0].area();
    double max_area = trgs[0].area();
    for(int i = 0; i < ntrgs; ++i)
    {
        min_area = (trgs[i].area() < min_area) ? trgs[i].area() : min_area;
        max_area = (trgs[i].area() > max_area) ? trgs[i].area() : max_area;
    }

    std::cout << "\tMin Area : " << min_area << std::endl;
    std::cout << "\tMax Area : " << max_area << std::endl;
}

/**
 * Smooth Mesh
**/
void Mesh::smooth(double smoothing_intensity,size_t niter) {
    std::vector< intSet > neighbors(npts);
    for(int i = 0; i < npts; ++i)
    {
        for(intSet::iterator it = links[i].begin(); it != links[i].end(); ++it)
        {
            Triangle nt = getTrg(*it);
            for(size_t k = 0; k < 3; ++k)
            {
                if (nt[k] != i) neighbors[i].insert(nt[k]);
            }
        }
    }
    Vect3* new_pts = new Vect3[npts];

    for(size_t n = 0; n < niter; ++n)
    {
        for(int p = 0; p < npts; ++p)
        {
            new_pts[p] = pts[p];
            for(intSet::iterator it = neighbors[p].begin(); it != neighbors[p].end(); ++it)
            {
                new_pts[p] = new_pts[p] + (smoothing_intensity * (pts[*it] - pts[p])) / neighbors[p].size();
            }
        }
        for(int p = 0; p < npts; ++p) pts[p] = new_pts[p];
    }
    delete[] new_pts;
    update_triangles(); // Updating triangles (areas + normals)
    recompute_normals(); // Updating normals
}

/**
 * Surface Gradient
**/
SparseMatrix Mesh::gradient() const {
    SparseMatrix A(3*ntrgs,npts); // nb edges x points
    // loop on triangles
    for (int t=0;t<ntrgs;t++) {
        const Triangle& trg = getTrg(t);
        Vect3 pts[3] = {getPt(trg[0]), getPt(trg[1]), getPt(trg[2])};
        for(int j=0;j<3;j++) {
            Vect3 grads = P1Vector(pts[0], pts[1], pts[2], j);
            for(int i=0;i<3;i++) {
                A(3*t+i,trg[j]) = grads(i);
            }
        }
    }
    return A;
}

Vector Mesh::areas() const {
    Vector areas(nbTrgs());
    for(int i = 0; i < nbTrgs(); ++i)
    {
        areas(i) = getTrg(i).getArea();
    }
    return areas;
}

void Mesh::update_triangles() {
    for(int i = 0; i < ntrgs; ++i)
    {
        trgs[i].normal() = pts[trgs[i][0]].normal( pts[trgs[i][1]] , pts[trgs[i][2]] );
        trgs[i].area() = trgs[i].normal().norm()/2.0;
    }
}

void Mesh::recompute_normals() {
    for(int p = 0; p < nbPts(); ++p)
    {
        Vect3 normal(0);
        for(intSet::iterator it = links[p].begin(); it != links[p].end(); ++it)
        {
            normal += trgs[*it].normal().normalize();
        }
        normal.normalize();
        normals[p] = normal;
    }
}

bool Mesh::selfIntersection() const {
    bool selfIntersects = false;
    for(int i = 0; i < ntrgs; ++i)
    {
        const Triangle& T1 = getTrg(i);
        for(int j = i+1; j < ntrgs; ++j)
        {
            const Triangle& T2 = getTrg(j);
            if(!T1.contains(T2.s1()) && !T1.contains(T2.s2()) && !T1.contains(T2.s3()))
            {
                if(triangle_intersection(*this,i,*this,j)) {
                    selfIntersects = true;
                }
            }
        }
    }
    return selfIntersects;
}

bool Mesh::intersection(const Mesh& m) const {
    bool intersects = false;
    for(int i = 0; i < ntrgs; ++i)
    {
        for(int j = 0; j < m.nbTrgs(); ++j)
        {
            intersects = intersects | triangle_intersection(*this,i,m,j);
        }
    }
    return intersects;
}
