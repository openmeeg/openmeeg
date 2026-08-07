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

#include "commandline.h"

using namespace OpenMEEG;

int main(int argc,char *argv[]) try {

    CommandLine cmd(argc,argv,"Add an array to a vtk file (attached to points or cells depending on the -c or -p options,\n automatically to one of them if dimension allows disambiguation).");
    const std::string& input_mesh  = cmd.option("-m",std::string(),"Input VTK/VTP polydata file");
          std::string  output_mesh = cmd.option("-o",std::string(),"Optional output filename for the vtk file");
    const std::string& data_name   = cmd.option("-n",std::string(),"Data field name");
          bool         point_type  = cmd.option("-p",false,"Matrix data represents point data");
          bool         cell_type   = cmd.option("-c",false,"Matrix data represents cell data");
    const std::string& data_file   = cmd.positional("data_file",std::string(),"Data file name");

    if (cmd.help_mode()) {
        return 0;
    }

    if (argc<6) {
        std::cerr << "Wrong number of arguments !" << std::endl;
        return 1;
    }

    if (input_mesh.empty()) {
        std:cerr << "A vtk input mesh must be given with -m !" << std::endl;
        return 1;
    }

    if (data_file.empty()) {
        std::cerr << "A data file must be given !" << std::endl;
        return 2;
    }

    // Read the vtk file (vtkPolyData).

    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(input_mesh.c_str());
    vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
    reader->Update();

    // Load the data.

    Matrix data(data_name);
    if (data.ncol()==0 || data.nlin()==0) {
        std::cerr << "You must give a correct matrix data file !" << std::endl;
        return 3;
    }

    if (output_mesh.empty())
        output_mesh = input_mesh;

    // Get the number of cells and points so we can know where to add our data

    const vtkIdType nbPoints = mesh->GetNumberOfPoints();
    const vtkIdType nbCells  = mesh->GetNumberOfCells();

    if (point_type && cell_type) {
        std::cerr << "-p and -c options are mutually exclusive. Please choose one or the other." << std::endl;
        return 4;
    }

    if (point_type) {
        if (data.nlin()!=static_cast<size_t>(nbPoints)) {
            std::cerr << "-p option is specified but the provided number of lines of the provided array does not match the number of points in the mesh (" << data.nlin() << " vs " << nbPoints << ")." << std::endl;
            return 5;
        }
    } else if (cell_type) {
        if (data.nlin()!=static_cast<size_t>(nbCells)) {
            std::cerr << "-c option is specified but the provided number of lines of the provided array does not match the number of cells in the mesh (" << data.nlin() << " vs " << nbCells << ")." << std::endl;
            return 5;
        }
    } else {
        if (nbPoints==nbCells) {
            std::cerr << "The provided vtk mesh has the same number of points and cells." << std::endl
                      << "It is thus impossible to infer to which type of entity the data has to be attached." << std::endl
                      << "Use the -p or -c options to select points or cells." << std::endl;
            return 6;
        }

        if (data.nlin()==static_cast<size_t>(nbPoints))
            point_type = true;
        else if (data.nlin()==static_cast<size_t>(nbCells))
            cell_type = true;
        else {
            std::cerr << "Something's wrong ! The input mesh has " << nbPoints << " points and " << nbCells << " cells, and the data matrix has "<< data.nlin() << " rows," << std::endl
                      << "which does not match any of these numbers ! Impossible attach such data to the mesh." << std::endl;
            return 7;
        }
    }

    // Add the data set   

    for (unsigned j=0; j<data.ncol(); ++j) {
        std::ostringstream array_name;
        array_name << data_name;

        // Build a vtkDataArray with our data

        vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
        for (unsigned i=0; i<data.nlin(); ++i)
            array->InsertNextValue(data(i,j));

        if (data.ncol()>1)
            array_name << '-' << std::setw(unsigned(log10(data.ncol())+1)) << std::setfill('0') << j;

        array->SetName(array_name.str().c_str());

        if (point_type)
            mesh->GetPointData()->AddArray(array);
        else if (cell_type)
            mesh->GetCellData()->AddArray(array);
    }

    // Systematically add indices arrays to points and cells.
    // If the arrays already exist (with the name "Indices"), this does nothing.

    vtkSmartPointer<vtkUnsignedIntArray> cell_indices = vtkSmartPointer<vtkUnsignedIntArray>::New();
    cell_indices->SetName("Indices");
    for (unsigned i=nbPoints; i<nbPoints+nbCells; ++i)
        cell_indices->InsertNextValue(i);
    mesh->GetCellData()->AddArray(cell_indices);

    vtkSmartPointer<vtkUnsignedIntArray> point_indices = vtkSmartPointer<vtkUnsignedIntArray>::New();
    point_indices->SetName("Indices");
    for (unsigned i=0; i<nbPoints; ++i)
        point_indices->InsertNextValue(i);
    mesh->GetPointData()->AddArray(point_indices);

    // Save output to the input file.

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(output_mesh.c_str());
    writer->SetInputData(mesh);
    writer->Write();

    std::cerr << "Saved into " << output_mesh << "." << std::endl;

    return 0;
} catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
}
