#ifndef H_fast_sparse_matrice
#define H_fast_sparse_matrice

#include <cstring>
#include <fstream>

#include "MatLibConfig.h"
#include "vecteur.h"
#include "sparse_matrice.h"

class fast_sparse_matrice : public genericMatrix
{
public:
    typedef sparse_matrice::idxType idxType;
    typedef sparse_matrice::valType valType;

protected:

    valType *tank;
    idxType *js;
    idxType *rowindex;
    idxType m_nlin;
    idxType m_ncol;

    inline void alloc(idxType nl, idxType nc, idxType nz);
    inline void destroy();

public:
    inline fast_sparse_matrice();
    inline fast_sparse_matrice(idxType n,idxType p, idxType sp);
    inline fast_sparse_matrice( const sparse_matrice &M);
    inline fast_sparse_matrice( const fast_sparse_matrice &M);
    inline ~fast_sparse_matrice() {destroy();}
    inline size_t nlin() const ;
    inline size_t ncol() const ;
    inline void saveTxt( const char *filename ) const;
    inline void saveBin( const char *filename ) const;
    inline void loadTxt( const char *filename );
    inline void loadBin( const char *filename );
    inline void write(std::ostream& f) const;
    inline void read(std::istream& f);

    inline double* getbuf() const;
    inline double operator()(size_t i,size_t j) const;
    inline double& operator()(size_t i,size_t j);
    inline vecteur operator * (const vecteur &v) const;
    inline void operator =( const fast_sparse_matrice &M);
    inline vecteur solve(const vecteur &rhs_vec) const;

    inline friend std::ostream& operator<<(std::ostream& f,const fast_sparse_matrice &M);
};

inline fast_sparse_matrice::fast_sparse_matrice()
{
    alloc(1,1,1);
}

inline fast_sparse_matrice::fast_sparse_matrice(idxType n,idxType p, idxType sp=1)
{
    alloc(n,p,sp);
}

inline void fast_sparse_matrice::operator =( const fast_sparse_matrice &M)
{
    destroy();
    alloc(M.m_nlin,M.m_ncol,M.rowindex[M.m_nlin]);
    memcpy(tank,M.tank,sizeof(valType)*M.rowindex[M.m_nlin]);
    memcpy(js,M.js,sizeof(idxType)*M.rowindex[M.m_nlin]);
    memcpy(rowindex,M.rowindex,sizeof(idxType)*(M.m_nlin+1));
}

inline fast_sparse_matrice::fast_sparse_matrice( const sparse_matrice &M)
{
    tank=new valType[M.getNz()];
    js=new idxType[M.getNz()];
    rowindex=new idxType[M.nlin()+1];
    m_nlin=(idxType)M.nlin();
    m_ncol=(idxType)M.ncol();

    // On remplit une structure plus rapide pour le calcul.
    sparse_matrice::cellType *cell;
    int cpt=0;
    for(idxType i=0;i<m_nlin;i++)
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
    rowindex[m_nlin]=cpt;

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
    ifs>>m_nlin>>m_ncol;
    ifs>>nz;
    alloc(m_nlin,m_ncol,nz);
    for(idxType i=0;i<nz;i++) ifs>>tank[i];
    for(idxType i=0;i<nz;i++) ifs>>js[i];
    for(idxType i=0;i<m_nlin+1;i++) ifs>>rowindex[i];

}
inline void fast_sparse_matrice::loadBin( const char *filename )
{
    std::ifstream ifs(filename,std::ios_base::binary);
    read(ifs);
    ifs.close();
}

inline void fast_sparse_matrice::write(std::ostream& f) const
{
    idxType nz=rowindex[m_nlin];
    f.write((const char*)&m_nlin,(std::streamsize)sizeof(idxType));
    f.write((const char*)&m_ncol,(std::streamsize)sizeof(idxType));
    f.write((const char*)&nz,(std::streamsize)sizeof(idxType));
    f.write((const char*)tank,(std::streamsize)(sizeof(valType)*nz));
    f.write((const char*)js,(std::streamsize)(sizeof(idxType)*nz));
    f.write((const char*)rowindex,(std::streamsize)(sizeof(idxType)*m_nlin));
}

inline void fast_sparse_matrice::read(std::istream& f)
{
    destroy();
    idxType nz;
    f.read((char*)&m_nlin,(std::streamsize)sizeof(idxType));
    f.read((char*)&m_ncol,(std::streamsize)sizeof(idxType));
    f.read((char*)&nz,(std::streamsize)sizeof(idxType));
    alloc(m_nlin,m_ncol,nz);
    f.read((char*)tank,(std::streamsize)(sizeof(valType)*nz));
    f.read((char*)js,(std::streamsize)(sizeof(idxType)*nz));
    f.read((char*)rowindex,(std::streamsize)(sizeof(idxType)*m_nlin));
}

inline std::ostream& operator<<(std::ostream& f,const fast_sparse_matrice &M) 
{
    for(fast_sparse_matrice::idxType i=0;i<M.m_nlin;i++)
    {
        for(fast_sparse_matrice::idxType j=M.rowindex[i];j<M.rowindex[i+1];j++)
        {
            f<<(long unsigned int)i<<"\t"<<(long unsigned int)M.js[j]<<"\t"<<M.tank[j];
        }
    }
    return f;
}

inline void fast_sparse_matrice::alloc(idxType nl, idxType nc, idxType nz)
{
    m_nlin=nl;
    m_ncol=nc;
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
    alloc(M.m_nlin,M.m_ncol,M.rowindex[M.m_nlin]);
    memcpy(tank,M.tank,sizeof(valType)*M.rowindex[M.m_nlin]);
    memcpy(js,M.js,sizeof(idxType)*M.rowindex[M.m_nlin]);
    memcpy(rowindex,M.rowindex,sizeof(idxType)*(M.m_nlin+1));
}

inline size_t fast_sparse_matrice::nlin() const {return (size_t)m_nlin;}

inline size_t fast_sparse_matrice::ncol() const {return (size_t)m_ncol;}

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
}

inline vecteur fast_sparse_matrice::operator * (const vecteur &v) const
{
    vecteur result(m_nlin); result.set(0);
    double *pt_result=&result(0);
    vecteur *_v=(vecteur *)&v;
    double *pt_vect=&(*_v)(0);

    for(idxType i=0;i<m_nlin;i++)
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
