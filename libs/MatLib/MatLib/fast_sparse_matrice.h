/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
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

    inline friend std::ostream& operator<<(std::ostream& f,const fast_sparse_matrice &M);

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

    inline double& operator[](size_t i) {return tank[i];};

    virtual std::ostream& operator>>(std::ostream& f) const
    {
        idxType nz = rowindex[m_nlin];
        f << m_nlin << " " << m_ncol << std::endl;
        f << nz << std::endl;
        for(idxType i=0;i<m_nlin;i++)
        {
            for(idxType j=rowindex[i];j<rowindex[i+1];j++)
            {
                f<<(long unsigned int)i<<"\t"<<(long unsigned int)js[j]<<"\t"<<tank[j]<<std::endl;
            }
        }
        return f;
    }
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
    idxType nz = rowindex[m_nlin];
    std::ofstream ofs(filename);
    ofs << m_nlin << " " << m_ncol << " " << nz << " ";
    for(idxType i=0;i<nz;i++) ofs << tank[i] << " ";
    for(idxType i=0;i<nz;i++) ofs << js[i] << " ";
    for(idxType i=0;i<m_nlin+1;i++) ofs << rowindex[i]  << " ";
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
    if(!ifs.is_open()) {
        std::cerr<<"Error Opening Matrix File "<<filename<<std::endl;
        exit(1);
    }
    ifs >> m_nlin >> m_ncol;
    ifs >> nz;
    alloc(m_nlin,m_ncol,nz);
    for(idxType i=0;i<nz;i++) ifs>>tank[i];
    for(idxType i=0;i<nz;i++) ifs>>js[i];
    for(idxType i=0;i<m_nlin+1;i++) ifs>>rowindex[i];
}

inline void fast_sparse_matrice::loadBin( const char *filename )
{
    std::ifstream ifs(filename,std::ios_base::binary);
	if(!ifs.is_open()) {
        std::cerr<<"Error Opening Matrix File "<<filename<<std::endl;
        exit(1);
    }
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

inline void fast_sparse_matrice::alloc(idxType nl, idxType nc, idxType nz)
{
    m_nlin=nl;
    m_ncol=nc;
    tank=new valType[nz];
    js=new idxType[nz];
    rowindex=new idxType[nl+1];
    rowindex[nl]=nz;
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
        double& total=pt_result[i];
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
