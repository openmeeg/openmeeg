#include <iostream>
#include <OMassert.H>

#include <mesh.h>

using namespace OpenMEEG;

int are_equal(const Vertex& v1,const Vertex& v2,double tol=1e-12) {
    return (v1-v2).norm2()<tol;
}

int are_equal(const Mesh& m1,const Mesh& m2,double tol=1e-12) {
    if ((m1.vertices().size()!=m2.vertices().size()) || (m1.triangles().size()!=m2.triangles().size()))
        return 0;

    for (auto vit1=m1.vertices().begin(),vit2=m2.vertices().begin();vit1!=m1.vertices().end();vit1++,vit2++)
        if (!are_equal(**vit1,**vit2,tol))
            return 0;

    for (Triangles::const_iterator tit1=m1.triangles().begin(),tit2=m2.triangles().begin();tit1!=m1.triangles().end();++tit1,++tit2)
        for (Triangle::const_iterator sit1=tit1->begin(),sit2=tit2->begin();sit1!=tit1->end();++sit1,++sit2)
            if (!are_equal(**sit1,**sit2))
                return 0;
    return 1;
}


int main(int argc,char** argv) {
	if (argc!=2) {
        std::cerr << "Wrong nb of parameters" << std::endl;
        exit(1);
    }

    std::cout << "Mesh : " << argv[1] << std::endl;

    Mesh mesh_orig;
    mesh_orig.load(argv[1]);

    Mesh mesh;

    // // TRI

    mesh_orig.save("tmp.tri");
    mesh.load("tmp.tri");

    om_error(are_equal(mesh, mesh_orig));

    // VTK
    mesh_orig.save("tmp.vtk");
#ifdef USE_VTK
    mesh.load("tmp.vtk");
    om_error(are_equal(mesh, mesh_orig));
#endif

#ifdef USE_GIFTI
    // GIFTI

    mesh_orig.save("tmp.gii");
    mesh.load("tmp.gii");
    om_error(are_equal(mesh, mesh_orig));
#endif

    // MESH

    mesh_orig.save("tmp.mesh");
    mesh.load("tmp.mesh");
    om_error(are_equal(mesh, mesh_orig));

    // BND && OFF that do not store normals

    mesh_orig.save("tmp.bnd");
    mesh_orig.save("tmp.off");
    mesh_orig.info();

    Mesh mesh1;
    Mesh mesh2;
    mesh1.load("tmp.bnd");
    mesh2.load("tmp.off");
    om_error(are_equal(mesh1, mesh2));

    return 0;
}
