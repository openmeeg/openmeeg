#include <iostream>
#include <OMassert.H>

#include <mesh.h>

using namespace OpenMEEG;

int are_equal(const Vertex& v1, const Vertex& v2, double tol = 1e-12) {
    return (v1 - v2).norm2() < tol;
}

int are_equal(const Mesh& m1, const Mesh& m2, double tol = 1e-12) {
    if ( m1.nb_vertices() != m2.nb_vertices() ) {
        return 0;
    }
    if ( m1.nb_triangles() != m2.nb_triangles() ) {
        return 0;
    }
    for ( Mesh::const_vertex_iterator vit1 = m1.vertex_begin(), vit2 = m2.vertex_begin(); vit1 != m1.vertex_end(); vit1++, vit2++)
    {
        if ( !are_equal(**vit1, **vit2, tol)) {
            return 0;
        }
    }
    for ( Mesh::const_iterator tit1 = m1.begin(), tit2 = m2.begin(); tit1 != m1.end(); ++tit1, ++tit2)
    {
        for ( Triangle::const_iterator sit1 = tit1->begin(), sit2 = tit2->begin(); sit1 != tit1->end(); ++sit1, ++sit2) {
            if ( !are_equal(**sit1, **sit2) ) {
                return 0;
            }
        }
    }
    return 1;
}


int main (int argc, char** argv)
{
	if ( argc != 2 ) 
    {
        std::cerr << "Wrong nb of parameters" << std::endl;
        exit(1);
    }

    std::cout << "Mesh : " << argv[1] << std::endl;

    Mesh mesh;
    mesh.load(argv[1]);

    //  TRI

    mesh.save("tmp.tri");
    Mesh trimesh;
    trimesh.load("tmp.tri");

    om_error(are_equal(mesh,trimesh));

    // VTK
    mesh.save("tmp.vtk");
#ifdef USE_VTK
    Mesh vtkmesh;
    vtkmesh.load("tmp.vtk");
    om_error(are_equal(mesh,vtkmesh));
#endif

    // GIFTI
#ifdef USE_GIFTI
    mesh.save("tmp.gii");
    Mesh giftimesh;
    giftimesh.load("tmp.gii");
    om_error(are_equal(mesh,giftimesh));
#endif

    // MESH

    mesh.save("tmp.mesh");
    Mesh mmesh;
    mmesh.load("tmp.mesh");
    om_error(are_equal(mesh,mmesh));

    // BND && OFF that do not store normals

    mesh.save("tmp.bnd");
    mesh.save("tmp.off");
    mesh.info();

    Mesh mesh1;
    mesh1.load("tmp.bnd");

    Mesh mesh2;
    mesh2.load("tmp.off");
    om_error(are_equal(mesh1, mesh2));

    return 0;
}
