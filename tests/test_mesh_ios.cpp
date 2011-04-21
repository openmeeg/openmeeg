#include <iostream>
#include <assert.h>

#include "mesh3.h"

using namespace OpenMEEG;

int are_equal(const Vect3& v1, const Vect3& v2, double tol=1e-12) {
    return (v1 - v2).norm2() < tol;
}

int are_equal(const Mesh& M1, const Mesh& M2, double tol=1e-12) {
    if (M1.nbPts() != M2.nbPts()) {
        return 0;
    }
    if (M1.nbTrgs() != M2.nbTrgs()) {
        return 0;
    }
    for(int i = 0; i < M1.nbPts(); ++i) {
        if (!are_equal(M1.point(i), M2.point(i), tol)) {
            return 0;
        }
    }
    for(int i = 0; i < M1.nbTrgs(); ++i) {
        if (M1.triangle(i) != M2.triangle(i)) {
            return 0;
        }
    }
    for(int i = 0; i < M1.nbPts(); ++i) {
        if (!are_equal(M1.normal(i), M2.normal(i), tol)) {
            return 0;
        }
    }
    return 1;
}


int main (int argc, char** argv)
{
	if(argc!=2) {
        std::cerr << "Wrong nb of parameters" << std::endl;
        exit(1);
    }

    std::cout << "Mesh : " << argv[1] << std::endl;

    Mesh mesh_orig;
    mesh_orig.load(argv[1]);

    Mesh mesh = mesh_orig;

    // // TRI
    mesh.save("tmp.tri");
    mesh.load("tmp.tri");

    assert(are_equal(mesh, mesh_orig));

    // VTK
    mesh.save("tmp.vtk");
#ifdef USE_VTK
    mesh.load("tmp.vtk");
    assert(are_equal(mesh, mesh_orig));
#endif

    // MESH
    mesh.save("tmp.mesh");
    mesh.load("tmp.mesh");
    assert(are_equal(mesh, mesh_orig));

    // BND && OFF that do not store normals
    mesh.save("tmp.bnd");
    mesh.save("tmp.off");

    Mesh mesh1;
    Mesh mesh2;
    mesh1.load("tmp.bnd");
    mesh2.load("tmp.off");
    assert(are_equal(mesh1, mesh2));

    return 0;
}
