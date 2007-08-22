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
#include <iostream>
#include <fstream>
#include <sstream>

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
    // store in streamFormat the format of the Mesh :
    char extension[128];
    getNameExtension(filename,extension);
    if(!strcmp(extension,"vtk") || !strcmp(extension,"VTK")) streamFormat = Mesh::VTK ;
    else if(!strcmp(extension,"tri") || !strcmp(extension,"TRI")) streamFormat = Mesh::TRI;
    else if(!strcmp(extension,"bnd") || !strcmp(extension,"BND")) streamFormat = Mesh::BND;
}

void Mesh::load(const char* name, bool verbose) {
    char extension[128];
    getNameExtension(name,extension);
    if(!strcmp(extension,"vtk") || !strcmp(extension,"VTK")) load_vtk(name);
    else if(!strcmp(extension,"mesh") || !strcmp(extension,"MESH")) load_Mesh(name);
    else if(!strcmp(extension,"tri") || !strcmp(extension,"TRI")) load_tri(name);
    else if(!strcmp(extension,"bnd") || !strcmp(extension,"BND")) load_bnd(name);
    else {
        cerr << "Load : Unknown file format" << endl;
        exit(1);
    }

    if(verbose)
    {
        std::cout<<"Mesh File: "<<name<<std::endl;
        std::cout<<"\t# points: "<<npts<<std::endl;
        std::cout<<"\t# triangles: "<<ntrgs<<std::endl;
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

    vtkDataArray *normalsData = vtkMesh->GetPointData()->GetNormals();
    assert(npts == normalsData->GetNumberOfTuples());
    assert(3 == normalsData->GetNumberOfComponents());
    
    for (int i = 0;i<npts;i++)
    {
        pts[i].x() = vtkMesh->GetPoint(i)[0];
        pts[i].y() = vtkMesh->GetPoint(i)[1];
        pts[i].z() = vtkMesh->GetPoint(i)[2];

        normals[i][0] = normalsData->GetTuple(i)[0];
        normals[i][1] = normalsData->GetTuple(i)[1];
        normals[i][2] = normalsData->GetTuple(i)[2];
    }
    ntrgs = vtkMesh->GetNumberOfCells();
    trgs = new Triangle[ntrgs];

    vtkIdList *l;
    for (int i = 0;i<ntrgs;i++) {
        if (vtkMesh->GetCellType(i) == VTK_TRIANGLE) {
            l = vtkMesh->GetCell(i)->GetPointIds();
            trgs[i][0] = l->GetId(0);
            trgs[i][1] = l->GetId(1);
            trgs[i][2] = l->GetId(2);
            trgs[i].normal() = pts[trgs[i][0]].normal( pts[trgs[i][1]] , pts[trgs[i][2]] );
            trgs[i].area() = trgs[i].normal().norme()/2.0;
        } else {
            std::cerr << "This is not a triangulation" << std::endl;
            exit(1);
        }
    }
}

void Mesh::load_vtk(std::istream &is) {
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
}

void Mesh::load_vtk(const char* name) {
    kill();

    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(name); // Specify file name of vtk data file to read
    if (!reader->IsFilePolyData()) {
        std::cerr << "This is not a valid vtk poly data file" << std::endl;
        reader->Delete();
        exit(1);
    }

    Mesh::getDataFromVTKReader(reader);
}

#endif

void Mesh::load_Mesh(std::istream &is) {
    std::cerr << "load_Mesh not inmplemented yet" << std::endl;
    exit(1);
}

void Mesh::load_Mesh(const char* name) {
    std::cerr << "load_Mesh not implemented yet (name: "<<name << " )." << std::endl;
    exit(1);
}

void Mesh::load_tri(std::istream &f) {
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
        trgs[i].normal() = pts[trgs[i][0]].normal( pts[trgs[i][1]] , pts[trgs[i][2]] );
        trgs[i].area() = trgs[i].normal().norme()/2.0;
    }
}

void Mesh::load_tri(const char* filename) {
    string s = filename;
    cout << "load_tri : " << filename << endl;
    std::ifstream f(filename);
    if(!f.is_open()) {
        cerr << "Error opening TRI file: " << filename << endl;
        exit(1);
    }
    Mesh::load_tri(f);
    f.close();
}

void Mesh::load_bnd(std::istream &f) {
    char myline[256];
    f.seekg( 0, std::ios_base::beg );
    f.getline(myline,256);
    f.getline(myline,256);

    // Vect3 normals;
    f.seekg( 0, std::ios_base::beg );
    f.getline(myline,256);
    f.getline(myline,256);

    // istringstream iss;
    string st;

    f.getline(myline,256);
    {
        istringstream iss(myline); iss >> st;
        assert(st == "NumberPositions=");
        iss >> npts;
    }

    f.getline(myline,256); // skip : "UnitPosition mm"

    f.getline(myline,256); {istringstream iss(myline); iss >> st;}
    assert(st == "Positions");

    pts = new Vect3[npts];
    links = new intSet[npts];

    for( size_t i = 0; i < size_t(npts); i += 1 )
    {
        f>>pts[i];
        pts[i] = pts[i];
    }
    f.getline(myline,256);
    f.getline(myline,256);
    {
        istringstream iss(myline); iss >> st;
        assert(st == "NumberPolygons=");
        iss >> ntrgs;
    }

    f.getline(myline,256);
    {
        istringstream iss(myline); iss >> st;
        assert(st == "TypePolygons=");
        iss >> st;
        assert(st == "3");
    }

    f.getline(myline,256); {istringstream iss(myline); iss >> st;}
    assert(st == "Polygons");

    trgs = new Triangle[ntrgs];
    for (int i=0;i<ntrgs;i++){
        f>>trgs[i];
        
        trgs[i].normal() = pts[trgs[i][0]].normal( pts[trgs[i][1]] , pts[trgs[i][2]] );
        trgs[i].area() = trgs[i].normal().norme()/2.0;
    }
}

void Mesh::load_bnd(const char* filename) {
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

void Mesh::save(const char* name) {
    char extension[128];
    getNameExtension(name,extension);
    if(!strcmp(extension,"vtk") || !strcmp(extension,"VTK")) save_vtk(name);
    else if(!strcmp(extension,"mesh") || !strcmp(extension,"MESH")) save_mesh(name);
    else if(!strcmp(extension,"bnd") || !strcmp(extension,"BND")) save_bnd(name);
    else if(!strcmp(extension,"tri") || !strcmp(extension,"TRI")) save_tri(name);
    else {
        cerr << "Save : Unknown file format" << endl;
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

void Mesh::save_mesh(const char* name) {
    std::cerr<<"save_mesh not implemented yet (name: "<<name << " )."<<std::endl;
    exit(1);
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
        newNormals[i + old_npts] = this->normal(i);
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
    return;
}

