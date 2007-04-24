#ifndef H_MESH
#define H_MESH

#include <istream>
#include <fstream>
#include <list>
#include "vect3.h"
#include "triangle.h"

#ifdef USE_VTK
#include <vtkPolyDataReader.h>
#endif

typedef std::list<int> intlist;

class mesh
{

private:
    int nb_pts, nb_trg;
    vect3 *pts;
    triangle *trg;
    intlist *links;


public:

    enum Filetype { VTK, TRI, TREE_D, BND } ;

    mesh ();
    mesh (int, int); // n_pts, n_trg
    mesh(int, vect3 *, int ,triangle* );
    mesh(const mesh& M);
    mesh& operator=( mesh& M);
    mesh& operator=(const mesh& M);
    ~mesh(){kill();}
    inline int nbr_pts() const { return nb_pts; }
    inline int nbr_trg() const { return nb_trg; }
    inline const vect3& vctr(int i) const { return pts[i]; }
    inline const triangle& trngl(int i) const { return trg[i]; }
    inline const intlist& trg_list(int i) const { return links[i]; }

    inline vect3& operator[](int i) {
        assert(i<nb_pts);
        return pts[i];
    }

    inline void set_pt(int i,const vect3& v) { pts[i]=v;}
    inline void set_trg(int i,const triangle& t) { trg[i]=t;}

    // each point knows each triangle it is a part of.
    void make_links();
    void compute_omega();
    vect3 center(int) const; // center of triangle i

    void load(const char* filename, bool verbose=true);

    void load_3d(std::istream &is);
    void load_3d(const char*);

#ifdef USE_VTK
    void load_vtk(std::istream &is);
    void load_vtk(const char*);
#else
    template <typename T>
    void load_vtk(T) {
        std::cerr << "You have to compile OpenMEEG with VTK to read VTK files" << std::endl;
        exit(1);
    }
#endif

    void load_tri(std::istream &);
    void load_tri(const char*);

    void load_bnd(std::istream &);
    void load_bnd(const char*);

    void load_mesh(std::istream &is);
    void load_mesh(const char*);

    void save(const char* filename);
    void save_3d(const char*);
    void save_vtk(const char*);
    void save_bnd(const char*);
    void save_tri(const char*);
    void save_mesh(const char*);

    inline int elem(int i, int* T ) const{ // gives each triangle a point belongs to
        int c=0;
        for(int k=0; k<nb_trg; k++){
            if( (i==trg[k].som1()) || (i==trg[k].som2()) || (i==trg[k].som3())){
                T[c] = k;
                c = c+1;
            }
        }
        return c;
    }

    void getFileFormat(const char* filename);

    inline friend void operator>>(std::istream &ifs,mesh &m){
        Filetype format = m.streamFormat;
        switch (format){
            case TRI:        m.load_tri(ifs); break;
            case TREE_D:     m.load_3d(ifs); break;
#ifdef USE_VTK            
            case VTK:        m.load_vtk(ifs); std::cout<<"Load vtk file format - OK "<< std::endl; break;
#endif 
            case BND:     m.load_bnd(ifs); /*std::cout << "mesh3.h - operator >> - switch 3D " << std::endl;*/ break;
            default: std::cout << "Format different de TRI ou 3d" << std::endl;
        }
    }

private:
    Filetype streamFormat;
    void kill();
    void copy(const mesh& M);
#ifdef USE_VTK
    void getDataFromVTKReader(vtkPolyDataReader* vtkMesh);
#endif    
};

#endif
