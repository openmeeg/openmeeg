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

#ifndef MESHDESCRIPTION_IO_H
#define MESHDESCRIPTION_IO_H

#include <iostream>
#include <fstream>
#include <vector>
#include <new>
#include <string>

#include <point.h>
#include <triangle.h>

namespace MeshDescription {

    void load_tri(std::istream &);
    void load_tri(const char*);

    void load_bnd(std::istream &);
    void load_bnd(const char*);

    void load_off(std::istream &);
    void load_off(const char*);

    void load_mesh(std::istream &is);
    void load_mesh(const char*);

#ifdef USE_VTK
    void load_vtk(std::istream &is);
    void load_vtk(const char*);
    void get_data_from_vtk_reader(vtkPolyDataReader* vtkMesh);
#else
    template <typename T>
        void load_vtk(T) {
            std::cerr << "You have to compile OpenMEEG with VTK to read VTK files" << std::endl;
            exit(1);
        }
#endif

#ifdef USE_GIFTI
    void load_gifti(std::istream &is);
    void load_gifti(const char*);
    void save_gifti(const char*);
    gifti_image* to_gifti_image();
    void from_gifti_image(gifti_image* gim);
#else
    template <typename T>
        void load_gifti(T) {
            std::cerr << "You have to compile OpenMEEG with GIFTI to read GIFTI files" << std::endl;
            exit(1);
        }
    template <typename T>
        void save_gifti(T) {
            std::cerr << "You have to compile OpenMEEG with GIFTI to save GIFTI files" << std::endl;
            exit(1);
        }
#endif

    void save_vtk(const char*) const;
    void save_bnd(const char*) const;
    void save_tri(const char*) const;
    void save_off(const char*) const;
    void save_mesh(const char*) const;

    Filetype streamFormat;

    void load(const char* filename) {

        std::string extension = getNameExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
        if (extension == std::string("vtk")) {
            load_vtk(filename);
        } else if (extension == std::string("tri")) {
            load_tri(filename);
        } else if (extension == std::string("bnd")) {
            load_bnd(filename);
        } else if (extension == std::string("mesh")) {
            load_mesh(filename);
        } else if (extension == std::string("off")) {
            load_off(filename);
        } else if (extension == std::string("gii")) {
            load_gifti(filename);
        } else {
            std::cerr << "IO: load: Unknown mesh file format for " << filename << std::endl;
            exit(1);
        }
    }


}
#endif  // ! MESHDESCRIPTION_IO_H
