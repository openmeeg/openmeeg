// Author: Jean-Christophe Lombardo
// Copyright (C) 2011-2017 - Jean-Christophe Lombardo, Inria

#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnsignedIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataWriter.h>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <matrix.h>


int usage(const std::string progname)
{
    std::cerr<<"Usage:"<<std::endl
             <<"\t"<<progname<<" -m <vtk polydata file> -n <data field name> [data file] (-o optional_output_filename) "<<std::endl
             <<""<<std::endl;
    return -1;
}

using namespace OpenMEEG;

int main(int argc, char *argv[])
{
    if (argc < 6) {
        std::cerr<<"Wrong number of arguments !"<<std::endl;
        return usage(argv[0]);
    }

    if ( (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") ) {
        return usage(argv[0]);
    }

    std::string meshFileName;
    std::string meshFileNameO = meshFileName;
    vtkSmartPointer<vtkPolyData> mesh;
    std::string dataName;
    Matrix data;

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
        } else if (std::string(argv[i]) == "-o") {
            // set an output filename
            meshFileNameO = argv[++i];
        } else {
            std::cerr<<"Opening data file <"<<argv[i]<<">"<<std::endl;
            data.load(argv[i]);
        }
    }

    if (! mesh.GetPointer()) {
        std::cerr<<"You must give a vtk mesh !"<<std::endl;
        return usage(argv[0]);
    }

    if ( data.ncol() * data.nlin() == 0) {
        std::cerr<<"You must give a correct matrix data file !"<<std::endl;
        return usage(argv[0]);
    }

    // get the number of cells and points so we can know where to add our data
    vtkIdType nbPoints = mesh->GetNumberOfPoints();
    vtkIdType nbCells = mesh->GetNumberOfCells();

    // Add the data set   
    for ( unsigned j = 0; j < data.ncol(); ++j) {
        std::ostringstream dataname;
        dataname << dataName;
        // Build a vtkDataArray with our data
        vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
        for ( unsigned i = 0; i < data.nlin(); ++i) {
            array->InsertNextValue(data(i, j));
        }
        if (data.ncol() > 1) {
            // dataname << "-" << std::setw(unsigned(log10(data.ncol())+1)) << std::setfill('0') << j; with 001 instead of 1
            dataname << "-" << j;
        }
        array->SetName(dataname.str().c_str());

        if (data.nlin() == (size_t)nbPoints) mesh->GetPointData()->AddArray(array);
        else if (data.nlin() == (size_t)nbCells) mesh->GetCellData()->AddArray(array);
        else {
            std::cerr << "Something's wrong ! got a mesh with " << nbPoints << " points and " << nbCells << " cells, and ("<< data.nlin() << ","<< data.ncol() <<") data. Don't know what to do with them !" << std::endl;
            return usage(argv[0]);
        }
    }

    // Add a indices dataset
    vtkSmartPointer<vtkUnsignedIntArray> cell_indices = vtkSmartPointer<vtkUnsignedIntArray>::New(); // indices
    vtkSmartPointer<vtkUnsignedIntArray> point_indices = vtkSmartPointer<vtkUnsignedIntArray>::New(); // indices
    cell_indices->SetName("Indices");
    point_indices->SetName("Indices");

    unsigned i = 0;
    for ( i = 0; i < nbPoints; ++i)
        point_indices->InsertNextValue(i);

    for (; i < nbPoints+nbCells; ++i)
        cell_indices->InsertNextValue(i);

    mesh->GetPointData()->AddArray(point_indices);
    mesh->GetCellData()->AddArray(cell_indices);

    // save output to the input file
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(meshFileNameO.c_str());
    #if VTK_MAJOR_VERSION >=6
    writer->SetInputData(mesh);
    #else
    writer->SetInput(mesh);
    #endif
    writer->Write();
    std::cerr << "saved into " << meshFileNameO << "." << std::endl;

    return 0;
}
