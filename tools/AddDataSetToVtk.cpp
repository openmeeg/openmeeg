// Author: Jean-Christophe Lombardo
// Copyright (C) 2011 - Jean-Christophe Lombardo, Inria

#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkGenericDataObjectWriter.h>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

int usage(const std::string progname)
{
    std::cerr<<"Usage:"<<std::endl
             <<"\t"<<progname<<" -m <vtk polydata file> -n <data field name> [data file]"<<std::endl
             <<""<<std::endl;
    return -1;
}

int main(int argc, char *argv[])
{
    if (argc != 6) {
        std::cerr<<"Wrong number of arguments !"<<std::endl;
        return usage(argv[0]);
    }

    if ( (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") ) {
        return usage(argv[0]);
    }


    std::string meshFileName;
    vtkSmartPointer<vtkPolyData> mesh;
    std::string dataName;
    float *data = 0l;
    unsigned long int nbData = 0;

    // parse command line arguments
    for (int i=1; i<argc; ++i) {
        if (std::string(argv[i]) == "-m") {
            // read a mesh (vtkPolyData)
            vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
            meshFileName = argv[++i];
            reader->SetFileName(meshFileName.c_str());
            mesh = reader->GetOutput();
            reader->Update();
        } else if (std::string(argv[i]) == "-n") {
            // read the new data field name
            dataName = argv[++i];
            std::cerr<<"DataName=<"<<dataName<<">"<<std::endl;
        } else {
            // read the data field
            // Not Yet implemented !
            std::cerr<<"Opening data file <"<<argv[i]<<">"<<std::endl;
            std::ifstream in(argv[i]);
            float val = 0;
            nbData = 0;
            // read the file a first time to count the data
            while (in.good()) {
                in>>val;
                if (in.good()) nbData++;
            }
            std::cerr<<"Got "<<nbData<<" data"<<std::endl;
            data = new float[nbData];
            
            in.clear(); // clear the error bit
            in.seekg(0, std::ios::beg); // rewind the file
            //in.close();
            //in.open();

            // read the file a second time to actually get the values
            unsigned long int idx=0;
            while (in.good()) {
              in>>val;
              data[idx++] = val;
            }
        }
    }

    if (! mesh.GetPointer()) {
        std::cerr<<"You must give a vtk mesh !"<<std::endl;
        return usage(argv[0]);
    }

    if (! data) {
        std::cerr<<"You must give a data file !"<<std::endl;
        return usage(argv[0]);
    }
 
    // Build a vtkDataArray with our data
    vtkSmartPointer<vtkFloatArray> array = vtkSmartPointer<vtkFloatArray>::New();
    array->SetArray(data, nbData, 0);
    array->SetName(dataName.c_str());

    // get the number of cells and points so we can know where to add our data
    vtkIdType nbPoints = mesh->GetNumberOfPoints();
    vtkIdType nbCells = mesh->GetNumberOfCells();

    if (nbData == nbPoints) mesh->GetPointData()->AddArray(array);
    else if (nbData == nbCells) mesh->GetCellData()->AddArray(array);
    else {
        std::cerr<<"Something's wrong ! got a mesh with "<<nbPoints<<" points and "<<nbCells<<" cells, and "<<nbData<<" data. Don't know what to do with them !"<<std::endl;
        return usage(argv[0]);
    }

    // save output to the input file
    vtkSmartPointer<vtkGenericDataObjectWriter> writer = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
    writer->SetFileName(meshFileName.c_str());
    writer->SetInput(mesh);
    writer->Write();

    return 0;
}
