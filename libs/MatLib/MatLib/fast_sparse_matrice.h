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

#include "MatLibConfig.h"
#include "vecteur_dcl.h"
#include "sparse_matrice_dcl.h"

class fast_sparse_matrice
{
public:

    inline friend std::ostream& operator<<(std::ostream& f,const fast_sparse_matrice &M);

protected:

    double *tank;
    size_t *js;
    size_t *rowindex;
    size_t m_nlin;
    size_t m_ncol;

    inline void alloc(size_t nl, size_t nc, size_t nz);
    inline void destroy();

public:
    inline fast_sparse_matrice();
    inline fast_sparse_matrice(size_t n,size_t p, size_t sp);
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

};

inline std::ostream& operator<<(std::ostream& f,const fast_sparse_matrice &M)
{
    size_t nz = M.rowindex[M.nlin()];
    f << M.nlin() << " " << M.ncol() << std::endl;
    f << nz << std::endl;
    for(size_t i=0;i<M.nlin();i++)
    {
        for(size_t j=M.rowindex[i];j<M.rowindex[i+1];j++)
        {
            f<<(long unsigned int)i<<"\t"<<(long unsigned int)M.js[j]<<"\t"<<M.tank[j]<<std::endl;
        }
    }
    return f;
}

inline fast_sparse_matrice::fast_sparse_matrice()
{
    alloc(1,1,1);
}

inline fast_sparse_matrice::fast_sparse_matrice(size_t n,size_t p, size_t sp=1)
{
    alloc(n,p,sp);
}

inline void fast_sparse_matrice::operator =( const fast_sparse_matrice &M)
{
    destroy();
    alloc(M.m_nlin,M.m_ncol,M.rowindex[M.m_nlin]);
    memcpy(tank,M.tank,sizeof(double)*M.rowindex[M.m_nlin]);
    memcpy(js,M.js,sizeof(size_t)*M.rowindex[M.m_nlin]);
    memcpy(rowindex,M.rowindex,sizeof(size_t)*(M.m_nlin+1));
}

inline fast_sparse_matrice::fast_sparse_matrice( const sparse_matrice &M)
{
    tank=new double[M.size()];
    js=new size_t[M.size()];
    rowindex=new size_t[M.nlin()+1];
    m_nlin=(size_t)M.nlin();
    m_ncol=(size_t)M.ncol();

    // we fill a data structure faster for computation
    sparse_matrice::const_iterator it;
    int cnt = 0;
    long current_line = -1;
    for( it = M.begin(); it != M.end(); ++it) {
        size_t i = it->first.first;
        size_t j = it->first.second;
        double val = it->second;
        tank[cnt] = val;
        js[cnt] = j;
        if(i != current_line) {
            for(size_t k = current_line+1; k <= i; ++k) {
                rowindex[k]=cnt;
            }
            current_line = i;
        }
        cnt++;
    }

    for(size_t k = current_line+1; k <= M.nlin(); ++k) {
        rowindex[k]=M.size();
    }

    this->saveTxt("totofspm.txt");
}

inline void fast_sparse_matrice::saveTxt( const char *filename ) const
{
    size_t nz = rowindex[m_nlin];
    std::ofstream ofs(filename);
    ofs << m_nlin << " " << m_ncol << " " << nz << " ";
    for(size_t i=0;i<nz;i++) ofs << tank[i] << " ";
    for(size_t i=0;i<nz;i++) ofs << js[i] << " ";
    for(size_t i=0;i<m_nlin+1;i++) ofs << rowindex[i]  << " ";
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
    size_t nz;
    std::ifstream ifs(filename);
    if(!ifs.is_open()) {
        std::cerr<<"Error Opening Matrix File "<<filename<<std::endl;
        exit(1);
    }
    ifs >> m_nlin >> m_ncol;
    ifs >> nz;
    alloc(m_nlin,m_ncol,nz);
    for(size_t i=0;i<nz;i++) ifs>>tank[i];
    for(size_t i=0;i<nz;i++) ifs>>js[i];
    for(size_t i=0;i<m_nlin+1;i++) ifs>>rowindex[i];
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
    size_t nz=rowindex[m_nlin];
    f.write((const char*)&m_nlin,(std::streamsize)sizeof(size_t));
    f.write((const char*)&m_ncol,(std::streamsize)sizeof(size_t));
    f.write((const char*)&nz,(std::streamsize)sizeof(size_t));
    f.write((const char*)tank,(std::streamsize)(sizeof(double)*nz));
    f.write((const char*)js,(std::streamsize)(sizeof(size_t)*nz));
    f.write((const char*)rowindex,(std::streamsize)(sizeof(size_t)*m_nlin));
}

inline void fast_sparse_matrice::read(std::istream& f)
{
    destroy();
    size_t nz;
    f.read((char*)&m_nlin,(std::streamsize)sizeof(size_t));
    f.read((char*)&m_ncol,(std::streamsize)sizeof(size_t));
    f.read((char*)&nz,(std::streamsize)sizeof(size_t));
    alloc(m_nlin,m_ncol,nz);
    f.read((char*)tank,(std::streamsize)(sizeof(double)*nz));
    f.read((char*)js,(std::streamsize)(sizeof(size_t)*nz));
    f.read((char*)rowindex,(std::streamsize)(sizeof(size_t)*m_nlin));
}

inline void fast_sparse_matrice::alloc(size_t nl, size_t nc, size_t nz)
{
    m_nlin=nl;
    m_ncol=nc;
    tank=new double[nz];
    js=new size_t[nz];
    rowindex=new size_t[nl+1];
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
    memcpy(tank,M.tank,sizeof(double)*M.rowindex[M.m_nlin]);
    memcpy(js,M.js,sizeof(size_t)*M.rowindex[M.m_nlin]);
    memcpy(rowindex,M.rowindex,sizeof(size_t)*(M.m_nlin+1));
}

inline size_t fast_sparse_matrice::nlin() const {return (size_t)m_nlin;}

inline size_t fast_sparse_matrice::ncol() const {return (size_t)m_ncol;}

inline double fast_sparse_matrice::operator()(size_t i,size_t j) const
{
    for(size_t k=rowindex[i];k<rowindex[i+1];k++)
    {
        if(js[k]<j) continue;
        else if(js[k]==j) return tank[k];
        else break;
    }

    return 0;
}

inline double& fast_sparse_matrice::operator()(size_t i,size_t j)
{
    for(size_t k=rowindex[i];k<rowindex[i+1];k++)
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

    for(size_t i=0;i<m_nlin;i++)
    {
        double& total=pt_result[i];
        for(size_t j=rowindex[i];j<rowindex[i+1];j++) {
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
