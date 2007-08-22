#ifndef H_fast_sparse_matrice
#define H_fast_sparse_matrice

#include <cstring>
#include <fstream>

#include "MatLibConfig.h"
#include "fast_sparse_matrice_dcl.h"
#include "vecteur.h"
#include "sparse_matrice.h"

inline fast_sparse_matrice::fast_sparse_matrice()
{
    alloc(1,1,1);
}

inline void fast_sparse_matrice::operator =( const fast_sparse_matrice &M)
{
    destroy();
    alloc(M.nlignes,M.ncolonnes,M.rowindex[M.nlignes]);
    memcpy(tank,M.tank,sizeof(valType)*M.rowindex[M.nlignes]);
    memcpy(js,M.js,sizeof(idxType)*M.rowindex[M.nlignes]);
    memcpy(rowindex,M.rowindex,sizeof(idxType)*(M.nlignes+1));
}

inline fast_sparse_matrice::fast_sparse_matrice( const sparse_matrice &M)
{
    tank=new valType[M.getNz()];
    js=new idxType[M.getNz()];
    rowindex=new idxType[M.nlin()+1];
    nlignes=(idxType)M.nlin();
    ncolonnes=(idxType)M.ncol();

    // On remplit une structure plus rapide pour le calcul.
    sparse_matrice::cellType *cell;
    int cpt=0;
    for(idxType i=0;i<nlignes;i++)
    {
        rowindex[i]=cpt;
        cell=((sparse_matrice*)&M)->getRowEntry()[i];
        if(cell!=0)
        {
            tank[cpt]=cell->val;
            js[cpt++]=cell->j;
            while(cell->right!=NULL)
            {
                cell=cell->right;
                tank[cpt]=cell->val;
                js[cpt++]=cell->j;
            }
        }
    }
    rowindex[nlignes]=cpt;

}

inline void fast_sparse_matrice::saveTxt( const char *filename ) const
{
    std::ofstream ofs(filename);
    ofs<<(*this);
    ofs.close();
}
inline void fast_sparse_matrice::saveBin( const char *filename ) const
{
    std::ofstream ofs(filename,std::ios_base::binary);
    write(ofs);
    ofs.close();
}
inline void fast_sparse_matrice::loadTxt( const char *filename )
{
    idxType nz;
    std::ifstream ifs(filename);
    ifs>>nlignes>>ncolonnes;
    ifs>>nz;
    alloc(nlignes,ncolonnes,nz);
    for(idxType i=0;i<nz;i++) ifs>>tank[i];
    for(idxType i=0;i<nz;i++) ifs>>js[i];
    for(idxType i=0;i<nlignes+1;i++) ifs>>rowindex[i];

}
inline void fast_sparse_matrice::loadBin( const char *filename )
{
    std::ifstream ifs(filename,std::ios_base::binary);
    read(ifs);
    ifs.close();
}

inline void fast_sparse_matrice::write(std::ostream& f) const
{
    idxType nz=rowindex[nlignes];
    f.write((const char*)&nlignes,(std::streamsize)sizeof(idxType));
    f.write((const char*)&ncolonnes,(std::streamsize)sizeof(idxType));
    f.write((const char*)&nz,(std::streamsize)sizeof(idxType));
    f.write((const char*)tank,(std::streamsize)(sizeof(valType)*nz));
    f.write((const char*)js,(std::streamsize)(sizeof(idxType)*nz));
    f.write((const char*)rowindex,(std::streamsize)(sizeof(idxType)*nlignes));
}

inline void fast_sparse_matrice::read(std::istream& f)
{
    destroy();
    idxType nz;
    f.read((char*)&nlignes,(std::streamsize)sizeof(idxType));
    f.read((char*)&ncolonnes,(std::streamsize)sizeof(idxType));
    f.read((char*)&nz,(std::streamsize)sizeof(idxType));
    alloc(nlignes,ncolonnes,nz);
    f.read((char*)tank,(std::streamsize)(sizeof(valType)*nz));
    f.read((char*)js,(std::streamsize)(sizeof(idxType)*nz));
    f.read((char*)rowindex,(std::streamsize)(sizeof(idxType)*nlignes));
}

inline std::ostream& operator<<(std::ostream& f,const fast_sparse_matrice &M) 
{
    for(fast_sparse_matrice::idxType i=0;i<M.nlignes;i++)
    {
        for(fast_sparse_matrice::idxType j=M.rowindex[i];j<M.rowindex[i+1];j++)
        {
            f<<(long unsigned int)i<<"\t"<<(long unsigned int)M.js[j]<<"\t"<<M.tank[j];
        }
    }
    return f;
}

//#if _MKL_ // TODO : a retablir quand le problem de link sera resolu
#if 0
vecteur fast_sparse_matrice::solve(const vecteur &rhs_vec) const
{
    vecteur result(nlignes);
    _DOUBLE_PRECISION_t *solValues=&result(0);
    _MKL_DSS_HANDLE_t handle;
    _INTEGER_t error;
    int opt = MKL_DSS_DEFAULTS;
    int non_sym = MKL_DSS_NON_SYMMETRIC;
    int type = MKL_DSS_INDEFINITE;

    const int un=1;
    int int_nlignes=(int)nlignes;
    int int_ncolonnes=(int)ncolonnes;
    int int_nz=(int)rowindex[nlignes];
    int *int_rowindex=new int[nlignes+1];
    for(int i=0;i<(int)nlignes+1;i++) {int_rowindex[i]=(int)rowindex[i]+1; std::cout<<int_rowindex[i]<<"  ";}; std::cout<<std::endl;
    int *int_js=new int[int_nz];
    for(int i=0;i<int_nz;i++) {int_js[i]=(int)js[i]+1; std::cout<<int_js[i]<<"  ";}; std::cout<<std::endl;
    std::cout<<1<<std::endl;

    error = dss_create(handle, opt ); //Initialize the solver
    if ( error != MKL_DSS_SUCCESS ) goto lbl_error; 

    std::cout<<2<<std::endl;
    error = dss_define_structure(
        handle, non_sym, int_rowindex, int_nlignes, int_ncolonnes,
        int_js, int_nz );
    if ( error != MKL_DSS_SUCCESS ) goto lbl_error;

    std::cout<<3<<std::endl;
    error = dss_reorder( handle, opt, 0); // Reorder the matrix
    if ( error != MKL_DSS_SUCCESS ) goto lbl_error;

    std::cout<<4<<std::endl;
    error = dss_factor_real( handle, type, tank ); // Factor the matrix
    if ( error != MKL_DSS_SUCCESS ) goto lbl_error;

    std::cout<<5<<std::endl;
    //&((((vecteur*)&rhs_vec)->operator()(0))
    error = dss_solve_real( handle, opt,&(((vecteur*)&rhs_vec)->operator()(0)), un, solValues ); // Get the solution
    if ( error != MKL_DSS_SUCCESS ) goto lbl_error;

    std::cout<<6<<std::endl;
    error = dss_delete( handle, opt ); // Deallocate solver storage
    if ( error != MKL_DSS_SUCCESS ) goto lbl_error;

    delete[] int_rowindex;
    delete[] int_js;
    return result;

    lbl_error:
    {
        delete[] int_rowindex;
        delete[] int_js;
        std::cerr<<"Error in fast_sparse_matrice : solve"<<std::endl;
        result.set(0.0);
        return result;
    }
}
#endif

inline void fast_sparse_matrice::alloc(idxType nl, idxType nc, idxType nz)
{
    nlignes=nl;
    ncolonnes=nc;
    tank=new valType[nz];
    js=new idxType[nz];
    rowindex=new idxType[nl+1];
}

inline void fast_sparse_matrice::destroy()
{
    delete[] tank;
    delete[] js;
    delete[] rowindex;
}

inline fast_sparse_matrice::fast_sparse_matrice( const fast_sparse_matrice &M)
{
    alloc(M.nlignes,M.ncolonnes,M.rowindex[M.nlignes]);
    memcpy(tank,M.tank,sizeof(valType)*M.rowindex[M.nlignes]);
    memcpy(js,M.js,sizeof(idxType)*M.rowindex[M.nlignes]);
    memcpy(rowindex,M.rowindex,sizeof(idxType)*(M.nlignes+1));
}

inline size_t fast_sparse_matrice::nlin() const {return (size_t)nlignes;}

inline size_t fast_sparse_matrice::ncol() const {return (size_t)ncolonnes;}

inline double fast_sparse_matrice::operator()(size_t i,size_t j) const
{
    for(idxType k=rowindex[i];k<rowindex[i+1];k++)
    {
        if(js[k]<j) continue;
        else if(js[k]==j) return tank[k];
        else break;
    }

    return 0;
}

inline double& fast_sparse_matrice::operator()(size_t i,size_t j)
{
    for(idxType k=rowindex[i];k<rowindex[i+1];k++)
    {
        if(js[k]<j) continue;
        else if(js[k]==j) return tank[k];
        else break;
    }

    std::cerr<<"fast_sparse_matrice : double& operator()(size_t i,size_t j) can't add element"<<std::endl;
    exit(1);
    //double *dummy=new double;
    //return *dummy; // devrait faire planter le programme
}

inline vecteur fast_sparse_matrice::operator * (const vecteur &v) const
{
    vecteur result(nlignes); result.set(0);
    double *pt_result=&result(0);
    vecteur *_v=(vecteur *)&v;
    double *pt_vect=&(*_v)(0);

    for(idxType i=0;i<nlignes;i++)
    {
        double& total=pt_result[i]; //deja à 0
        for(idxType j=rowindex[i];j<rowindex[i+1];j++) {
            total+=tank[j]*pt_vect[js[j]];
        }
    }
    return result;
}

inline double* fast_sparse_matrice::getbuf() const
{
    return tank;
}
#endif
