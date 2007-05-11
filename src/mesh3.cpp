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
#include <vtkCleanPolyData.h>
#endif

using namespace std;

mesh::mesh()
{
    nb_pts=0;
}

mesh::mesh(int a, int b)
{
    nb_pts=a;  nb_trg=b;
    pts = new vect3[nb_pts];
    trg = new triangle[nb_trg];
    links = new intlist[nb_pts];
}

mesh::mesh ( int a, vect3 *b, int c,triangle* d)
{
    nb_pts=a;
    pts=b;
    nb_trg=c;
    trg=d;
    links = new intlist[nb_pts];
}

mesh::mesh(const mesh& M){
    *this = M;
}

mesh& mesh::operator=(const mesh& M){
    if (this!=&M) copy(M);
    return *this;

}

//mesh& mesh::operator=(mesh& M){
//    if (this!=&M) copy(M);
//    return *this;
//}

void mesh::copy(const mesh& M){
    nb_pts = M.nb_pts;
    if (nb_pts!=0) {
        nb_trg = M.nb_trg;
        pts = new vect3[nb_pts];
        trg = new triangle[nb_trg];
        links = new intlist[nb_pts];
        for(int i=0;i<nb_pts;i++){
            pts[i] = M.pts[i];
            links[i] = M.links[i];
        }
        for(int i=0;i<nb_trg;i++)
            trg[i] = M.trg[i];
    }
}


void mesh::kill()
{
    if (nb_pts!=0) {
        delete []pts;
        delete []trg;
        delete []links;

        nb_pts=0;
        nb_trg=0;
    }
}

void mesh::make_links()
{
    int i;
    for (i=0;i<nb_pts;i++) {
        links[i].clear();
    }

    for (i=0;i<nb_trg;i++) {
        links[trg[i].som1()].push_back(i);
        links[trg[i].som2()].push_back(i);
        links[trg[i].som3()].push_back(i);
    }
}


void mesh::getFileFormat(const char* filename){
    //store in streamFormat the format of the mesh :
    char extension[128];
    getNameExtension(filename,extension);
    if(!strcmp(extension,"vtk") || !strcmp(extension,"VTK")) streamFormat=mesh::VTK ;
    else if(!strcmp(extension,"3d") || !strcmp(extension,"3D")) streamFormat=mesh::TREE_D;
    else if(!strcmp(extension,"tri") || !strcmp(extension,"TRI")) streamFormat=mesh::TRI;
    else if(!strcmp(extension,"bnd") || !strcmp(extension,"BND")) streamFormat=mesh::BND;
}


void mesh::load(const char* name, bool verbose) {
    char extension[128];
    getNameExtension(name,extension);
    if(!strcmp(extension,"vtk") || !strcmp(extension,"VTK")) load_vtk(name);
    else if(!strcmp(extension,"3d") || !strcmp(extension,"3D")) load_3d(name);
    else if(!strcmp(extension,"mesh") || !strcmp(extension,"MESH")) load_mesh(name);
    else if(!strcmp(extension,"tri") || !strcmp(extension,"TRI")) load_tri(name);
    else if(!strcmp(extension,"bnd") || !strcmp(extension,"BND")) load_bnd(name);
    else {
        cerr << "Load : Unknown file format" << endl;
        exit(1);
    }
    
    if(verbose)
    {
        std::cout<<"Mesh File: "<<name<<std::endl;
        std::cout<<"\t# points: "<<nb_pts<<std::endl;
        std::cout<<"\t# triangles: "<<nb_trg<<std::endl;
    }
}

#ifdef USE_VTK
void mesh::getDataFromVTKReader(vtkPolyDataReader* reader){   //private
    
    reader->Update();
    vtkPolyData *vtkMesh = reader->GetOutput();

    nb_pts = vtkMesh->GetNumberOfPoints();

    pts=new vect3[nb_pts]; // points array
    links=new intlist[nb_pts];

    for (int i=0;i<nb_pts;i++)
    {
        pts[i]._x() = vtkMesh->GetPoint(i)[0];
        pts[i]._y() = vtkMesh->GetPoint(i)[1];
        pts[i]._z() = vtkMesh->GetPoint(i)[2];
    }
    nb_trg = vtkMesh->GetNumberOfCells();
    trg = new triangle[nb_trg];

    vtkIdList *l;
    for (int i=0;i<nb_trg;i++) {
        if (vtkMesh->GetCellType(i) == VTK_TRIANGLE) {
            l = vtkMesh->GetCell(i)->GetPointIds();
            trg[i].som1() = l->GetId(0);
            trg[i].som2() = l->GetId(1);
            trg[i].som3() = l->GetId(2);
            trg[i].set_normale(pts[trg[i].som1()].normale(pts[trg[i].som2()],pts[trg[i].som3()]));
            trg[i].setAire(trg[i].normale().norme()/2.0);
        } else {
            std::cerr << "This is not a triangulation" << std::endl;
            exit(1);
        }
    }

}

void mesh::load_vtk(std::istream &is){
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
    
    mesh::getDataFromVTKReader(reader);  

    delete[] buffer;
    reader->Delete();
}



void mesh::load_vtk(const char* name) {
    kill();

    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(name); // Specify file name of vtk data file to read
    if (!reader->IsFilePolyData()) {
        std::cerr << "This is not a valid vtk poly data file" << std::endl;
        reader->Delete();
        exit(1);
    }
    
    mesh::getDataFromVTKReader(reader);
}

#endif



void mesh::load_3d(std::istream &in){
    kill();
    in>>nb_pts;
    pts=new vect3[nb_pts];
    links=new intlist[nb_pts];

    for(int i=0;i<nb_pts;i++){  // reading points
        in>>pts[i]._x();
        in>>pts[i]._y();
        in>>pts[i]._z();
    }

    in>>nb_trg;
    trg=new triangle[nb_trg]; // reading triangles
    for(int i=0;i<nb_trg;i++){
        in>>trg[i].som1();
        in>>trg[i].som2();
        in>>trg[i].som3();

        trg[i].set_normale(pts[trg[i].som1()].normale(pts[trg[i].som2()],pts[trg[i].som3()]));
        trg[i].setAire(trg[i].normale().norme()/2.0);
    }
}

void mesh::load_3d(const char *name){
    ifstream in(name,ios::in);
    if(!in.is_open()) {
        cerr<<"Error Reading File : " << name << std::endl;
        exit(1);
    }
    mesh::load_3d(in);
    in.close();
}

void mesh::load_mesh(std::istream &is){
    std::cerr << "load_mesh not inmplemented yet" << std::endl;
    exit(1);
}
void mesh::load_mesh(const char* name) {
    std::cerr << "load_mesh not implemented yet (name: "<<name << " )." << std::endl;
    exit(1);
}

void mesh::load_tri(std::istream &f){
    char ch; f>>ch;
    f>>nb_pts;

    char myline[256];
    f.seekg( 0, std::ios_base::beg );
    f.getline(myline,256);
    f.getline(myline,256);
    int nread=0;

    float r;
    nread = sscanf(myline,"%f %f %f %f %f %f",&r,&r,&r,&r,&r,&r);

    bool withNormals = (nread==6);
    vect3 normals;
    f.seekg( 0, std::ios_base::beg );
    f.getline(myline,256);

    pts=new vect3[nb_pts];
    links=new intlist[nb_pts];
    for (int i=0;i<nb_pts;i++){
        f>>pts[i];
        pts[i]=pts[i];
        if(withNormals) f >> normals; // Saute la normale (format Nicolas)
    }
    f>>ch;
    f>>nb_trg; f>>nb_trg; f>>nb_trg; // Ce numÈro est rÈpÈtÈ trois fois.
    trg=new triangle[nb_trg];
    for (int i=0;i<nb_trg;i++){
        f>>trg[i];
        trg[i].set_normale(pts[trg[i].som1()].normale(pts[trg[i].som2()],pts[trg[i].som3()]));
        trg[i].setAire(trg[i].normale().norme()/2.0);
    }
}

void mesh::load_tri(const char* filename) {
    string s=filename;
    cout << "load_tri : " << filename << endl;
    std::ifstream f(filename);
    if(!f.is_open()) {
        cerr << "Error opening TRI file: " << filename << endl;
        exit(1);
    }
    mesh::load_tri(f);
    f.close();
}
void mesh::load_bnd(std::istream &f){
    char myline[256];
    f.seekg( 0, std::ios_base::beg );
    f.getline(myline,256);
    f.getline(myline,256);

    // vect3 normals;
    f.seekg( 0, std::ios_base::beg );
    f.getline(myline,256);
    f.getline(myline,256);

    // istringstream iss;
    string st;

    f.getline(myline,256);
    {
        istringstream iss(myline); iss >> st;
        assert(st == "NumberPositions=");
        iss >> nb_pts;
    }

    f.getline(myline,256); // skip : "UnitPosition mm"

    f.getline(myline,256); {istringstream iss(myline); iss >> st;}
    assert(st == "Positions");

    pts = new vect3[nb_pts];
    links = new intlist[nb_pts];

    for( size_t i = 0; i < size_t(nb_pts); i += 1 )
    {
        f>>pts[i];
        pts[i]=pts[i];
    }
    f.getline(myline,256);
    f.getline(myline,256);
    {
        istringstream iss(myline); iss >> st;
        assert(st == "NumberPolygons=");
        iss >> nb_trg;
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

    trg=new triangle[nb_trg];
    for (int i=0;i<nb_trg;i++){
        f>>trg[i];
        trg[i].set_normale(pts[trg[i].som1()].normale(pts[trg[i].som2()],pts[trg[i].som3()]),true); // set normale with normalization
        // trg[i].set_normale(pts[trg[i].som1()].normale(pts[trg[i].som2()],pts[trg[i].som3()]),false); // set normale with normalization
        trg[i].setAire(trg[i].normale().norme()/2.0);
    }
}

void mesh::load_bnd(const char* filename) {
    string s=filename;
    cout << "load_bnd : " << filename << endl;
    std::ifstream f(filename);
    if(!f.is_open()) {
        cerr << "Error opening BND file: " << filename << endl;
        exit(1);
    }
    mesh::load_bnd(f);
    f.close();
}

void mesh::save(const char* name) {
    char extension[128];
    getNameExtension(name,extension);
    if(!strcmp(extension,"vtk") || !strcmp(extension,"VTK")) save_vtk(name);
    else if(!strcmp(extension,"3d") || !strcmp(extension,"3D")) save_3d(name);
    else if(!strcmp(extension,"mesh") || !strcmp(extension,"MESH")) save_mesh(name);
    else if(!strcmp(extension,"bnd") || !strcmp(extension,"BND")) save_bnd(name);
    else if(!strcmp(extension,"tri") || !strcmp(extension,"TRI")) save_tri(name);
    else {
        cerr << "Save : Unknown file format" << endl;
        exit(1);
    }
}

void mesh::save_3d(const char* name) {
    std::ofstream f(name);
    f<<nb_pts<<std::endl;
    for (int i=0;i<nb_pts;i++)
        f<<pts[i]<<std::endl;

    f<<nb_trg<<std::endl;
    for (int i=0;i<nb_trg;i++)
        f<<trg[i]<<std::endl;

    for (int i=0;i<nb_pts;i++)
        f<<trg[i].normale()<<std::endl;
    f.close();
}

void mesh::save_vtk ( const char* filename) {
    std::ofstream os(filename);
    os<<"# vtk DataFile Version 2.0"<<std::endl;
    os<<"File "<<filename<<" generated by OpenMEEG"<<std::endl;
    os<<"ASCII"<<std::endl;
    os<<"DATASET POLYDATA"<<std::endl;
    os<<"POINTS "<<nb_pts<<" float"<<std::endl;
    for(int i=0;i<nb_pts;i++) os<<pts[i]<<std::endl;
    os<<"POLYGONS "<<nb_trg<<" "<<nb_trg*4<<std::endl;
    for(int i=0;i<nb_trg;i++) os<<3<<" "<<trg[i]<<std::endl;
    os.close();
}

void mesh::save_bnd ( const char* filename) {
    std::ofstream os(filename);
    os << "# Bnd mesh file generated by OpenMeeg" << std::endl;
    os << "Type= Unknown" << std::endl;
    os << "NumberPositions= " << nb_pts << std::endl;
    os << "UnitPosition\tmm" << std::endl;
    os << "Positions" << std::endl;
    for(int i=0;i<nb_pts;i++) os<<pts[i]<<std::endl;
    os << "NumberPolygons= " << nb_trg <<std::endl;
    os << "TypePolygons=\t3" << std::endl;
    os << "Polygons" << std::endl;
    for(int i=0;i<nb_trg;i++) os << trg[i] << std::endl;
    os.close();
}

void mesh::save_tri ( const char* filename) {
    std::ofstream os(filename);
    os<<"- "<<this->nb_pts<<std::endl;
    for(int i=0;i<nb_pts;i++) {
        os<<this->pts[i]<<" ";
        os<<this->trg[i].normale()<<std::endl;
    }
    os<<"- "<<nb_trg<<" "<<nb_trg<<" "<<nb_trg<<std::endl;
    for(int i=0;i<nb_trg;i++) os<<this->trg[i]<<std::endl;
    os.close();
}

void mesh::save_mesh(const char* name) {
    std::cerr<<"save_mesh not implemented yet (name: "<<name << " )."<<std::endl;
    exit(1);
}

