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

#ifndef H_SPARSE_MATRICE
#define H_SPARSE_MATRICE

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <fstream>

#include "vecteur.h"

enum sparsity {row,column};
typedef sparsity sparsity;

template<class T,class I> struct MyCell
{
    T val;
    I i;
    I j;
    MyCell<T,I> *up;
    MyCell<T,I> *down;
    MyCell<T,I> *left;
    MyCell<T,I> *right;
};

template<class T,class I> class TsparseMatrixIterator;
template<class T,class I> class TsparseMatrix;

typedef TsparseMatrix<double,size_t> sparse_matrice;
typedef TsparseMatrixIterator<double,size_t> sparse_matrice_iterator;

template<class T,class I>
class TsparseMatrix
{
public:
    typedef T valType;
    typedef I idxType;
    typedef MyCell<T,I> cellType;

    friend class TsparseMatrixIterator<T,I>;

protected:
    // On a 2 tableaux d'indices pour parcourir les matrices
    cellType **_RowEntry;
    cellType **_ColEntry;
    bool _tr;

    sparsity _sp;

    // Nombre de lignes de la matrice
    idxType _n;
    // Nombre de colonnes de la matrice
    idxType _p;

    // Nombre d'elements non nuls dans la matrice
    idxType _nz;

    // tableaux des lignes contenant un element non nul
    idxType *_nzrows;
    idxType _nnzrows;

    // tableaux des colonnes contenant un element non nul
    idxType *_nzcols;
    idxType _nnzcols;

    inline void destroy()
    {
        for(idxType i=0;i<_n;i++)
        {
            cellType *pt,*ptbis;
            pt=_RowEntry[i];
            while(pt!=0)
            {
                ptbis=pt->right;
                delete pt;
                pt=ptbis;
            }
        }

        delete[] _RowEntry;
        delete[] _ColEntry;
        delete[] _nzrows;
        delete[] _nzcols;
    }

    inline void init()
    {
        _tr=false;
        _RowEntry=0;
        _ColEntry=0;
        _n=0;
        _p=0;
        _nz=0;
        _nnzrows=0;
        _nnzcols=0;
        _sp=column;
    }

    inline void init (idxType n,idxType p, sparsity sp=column)
    {
        _nz=0;
        _tr=false;

        _n=n;
        _p=p;
        _sp=sp;
        _RowEntry=new cellType*[n]; for(idxType k=0;k<_n;k++) _RowEntry[k]=idxType(0);
        _ColEntry=new cellType*[p]; for(idxType k=0;k<_p;k++) _ColEntry[k]=idxType(0);
        _nzrows=new idxType[n]; for(idxType i=0;i<_n;i++) _nzrows[i]=idxType(-1);
        _nzcols=new idxType[p]; for(idxType i=0;i<_p;i++) _nzcols[i]=idxType(-1);
        _nnzrows=n;
        _nnzcols=p;
    }

public:

    inline void refreshNZ()
    {
        _nnzrows=0;
        _nnzcols=0;
        idxType idx;
        idx=0;
        for(idxType i=0;i<_n;i++) if(_RowEntry[i]!=0) {_nnzrows++; _nzrows[idx]=i; idx++;}
        idx=0;
        for(idxType i=0;i<_p;i++) if(_ColEntry[i]!=0) {_nnzcols++; _nzcols[idx]=i; idx++;}
    }

    inline size_t nlin() const {return (size_t)_n;}
    inline size_t ncol() const {return (size_t)_p;}
    inline cellType **getColEntry() const {return _ColEntry;}
    inline cellType **getRowEntry() const {return _RowEntry;}
    inline idxType getNz() const {return _nz;}
    inline idxType getNnzrows() const { return _nnzrows;}
    inline idxType getNnzcols() const { return _nnzcols;}

    inline TsparseMatrix()
    {
        init();
    }

    inline TsparseMatrix(idxType n,idxType p, sparsity sp=column)
    {
        init(n,p,sp);
    }

    inline ~TsparseMatrix()
    {
        destroy();
    }

    inline TsparseMatrix ( const TsparseMatrix<T,I>& M )
    {
        _n=M._n;
        _p=M._p;
        //nz=M.nz;
        _nz=0;
        _sp=M._sp;
        _tr=M._tr;
        _RowEntry=new cellType*[_n]; for(idxType i=0;i<_n;i++) _RowEntry[i]=0;
        _ColEntry=new cellType*[_p]; for(idxType j=0;j<_p;j++) _ColEntry[j]=0;
        _nzrows=new idxType[_n]; for(idxType i=0;i<_n;i++) _nzrows[i]=M._nzrows[i];
        _nzcols=new idxType[_p]; for(idxType i=0;i<_p;i++) _nzcols[i]=M._nzcols[i];
        _nnzrows=M._nnzrows;
        _nnzcols=M._nnzcols;

        for(idxType i=0;i<_n;i++)
        {
            cellType *pt;
            pt = M._RowEntry[i];
            while(pt!=0)
            {
                (*this)(pt->i,pt->j) = M(pt->i,pt->j);
                pt=pt->right;
            }
        }
    }

    virtual std::ostream& operator>>(std::ostream& f) const
    {
        f << _n << " " << _p << std::endl;
        TsparseMatrixIterator<T,I> it(*this);
        for(it.begin(); !it.end(); it.next()) {
            f << (unsigned long int) it.current()->i << "\t";
            f << (unsigned long int) it.current()->j << "\t" << it.current()->val << std::endl;
        }
        return f;
    }

    inline cellType* first() const // get first non empty cell
    {
        cellType* cell;
        for (idxType i=0;i<_p;i++) {
            cell=_ColEntry[i];
            if (cell!=0) {
                return cell;
            }
        }
        return NULL;
    }

    vecteur operator*( const vecteur &x) const;

    inline TsparseMatrix<T,I> operator*(  TsparseMatrix<T,I>& opd )
    {

        refreshNZ();
        opd.refreshNZ();
        if(_p!=opd.nlin()) {std::cerr<<"Erreur Produit Matrice"; exit(1);}

        TsparseMatrix<T,I> result(_n,opd.ncol());

        for(idxType i=0;i<_nnzrows;i++)
            for(idxType j=0;j<opd._nnzcols;j++)
        {
            cellType* startCol=(opd._ColEntry)[opd._nzcols[j]];
            cellType* startRow=_RowEntry[_nzrows[i]];
            if(startCol!=0 && startRow!=0)
            {
                T total=0;
                while(startRow!=0)
                {
                    while( startCol->i < startRow->j && startCol->down!=0) startCol=startCol->down;
                        //if(startCol==0) goto finIJ;
                    if(startCol->i == startRow->j) total+=startCol->val*startRow->val;
                    startRow=startRow->right;
                }
                if (total!=0) (result)(_nzrows[i],opd._nzcols[j])=total;
            }
        }

        result.refreshNZ();
        return result;
    }

    inline void saveBin( const char * filename) const
    {
        FILE *file;
        file=fopen(filename,"wb");
	    if (file==NULL) {std::cerr<<"Error opening file: " << filename << std::endl; exit(1);}

        // Start by writting matrix dimensions
        idxType nn=this->_n;
        idxType pp=this->_p;
        fwrite(&nn,sizeof(idxType),1,file);
        fwrite(&pp,sizeof(idxType),1,file);

        // write each element
        TsparseMatrixIterator<T,I> it(*this);
        for(it.begin(); !it.end(); it.next()) {
            fwrite(&(it.current()->i),sizeof(idxType),1,file);
            fwrite(&(it.current()->j),sizeof(idxType),1,file);
            fwrite(&(it.current()->val),sizeof(valType),1,file);
        }

        fclose(file);
    }

    inline void saveTxt( const char * filename ) const
    {
        std::ofstream of(filename);
        *this >> of;
    }

    inline void saveMat( const char *filename ) const
    {
    #ifdef USE_MATIO
        mat_t* mat = Mat_Open(filename,MAT_ACC_RDWR);
        if (mat) {
            valType* t = (double*)malloc(sizeof(double)*_nz);
            mat_int32_t  *ir = (mat_int32_t*)malloc(sizeof(mat_int32_t)*_nz);
            mat_int32_t  *jc = (mat_int32_t*)malloc(sizeof(mat_int32_t)*(_p+1));

            jc[_p] = _nz;

            cellType *ce;

            size_t counter = 0;

            for(size_t j=0;j<_p;j++)
            {
                jc[j] = counter;
                if(_ColEntry[j] != NULL) // check if the column is not empty
                {
                    ce = _ColEntry[j];
                    do {
                        ir[counter] = ce->i;
                        t[counter] = ce->val;
                        ce = ce->down;
                        counter++;
                    } while (ce != 0);
                }
            }

            int dims[2] = { _n, _p };
            matvar_t *matvar;
            sparse_t  sparse = {0,};
            sparse.nzmax = _nz;
            sparse.nir   = _nz;
            sparse.ir    = ir;
            sparse.njc   = _p+1;
            sparse.jc    = jc;
            sparse.ndata = _nz;
            sparse.data  = t;

            matvar = Mat_VarCreate("matrix",MAT_C_SPARSE,MAT_T_DOUBLE,2,dims,&sparse,MEM_CONSERVE);
            Mat_VarWrite(mat,matvar,COMPRESSION_ZLIB);
            Mat_VarFree(matvar);
            Mat_Close(mat);
        }
    #else
        std::cerr << "You have to compile OpenMEEG with MATIO to save matlab files" << std::endl;
    #endif
    }

    inline void loadTxt( const char * filename ) // Les colonnes sont ecrites les unes à la suite des autres
    {
        load(filename,'t');
    }

    inline int readDim(FILE *file, idxType& nn, idxType& pp, const char c) {
        int nread = 0;
        if (c=='t') {
            int nnn;
            int ppp;
            nread += fscanf(file,"%d",&nnn);
            nread += fscanf(file,"%d",&ppp);
            nn = nnn;
            pp = ppp;
        } else if (c=='b') {
            nread += fread(&nn,sizeof(idxType),1,file);
            nread += fread(&pp,sizeof(idxType),1,file);
        } else {
            std::cerr << "Unknown format" << std::endl;
        }
        return nread;
    }

    inline int readCell(FILE *file, idxType& ii, idxType& jj, valType& val, const char c) {
        int nread = 0;
        if (c=='t') {
            int iii,jjj;
            nread += fscanf(file,"%d",&iii);
            nread += fscanf(file,"%d",&jjj);
            nread += fscanf(file,"%lf",&val);
            ii = iii;
            jj = jjj;
        } else if (c=='b') {
            nread += fread(&ii,sizeof(idxType),1,file);
            nread += fread(&jj,sizeof(idxType),1,file);
            nread += fread(&val,sizeof(valType),1,file);
        } else {
            std::cerr << "Unknown format" << std::endl;
        }
        return nread;
    }

    inline void loadBin( const char * filename) // Columns are written next to one another
    {
        load(filename,'b');
    }

    inline void loadMat( const char *filename ) throw(std::string)
    {
    #ifdef USE_MATIO
        mat_t* mat = Mat_Open(filename,MAT_ACC_RDONLY);
        if (mat) {
            matvar_t* matvar = Mat_VarReadNext(mat);

            while( matvar!=NULL && (matvar->rank!=2 || matvar->data_type!=MAT_T_DOUBLE || matvar->class_type!=MAT_C_SPARSE) )
                matvar = Mat_VarReadNext(mat);
            if (matvar==NULL)
                std::cerr << "There is no 2D sparse double matrix in file : " << filename << std::endl;

            std::cout << "loadMat : Using variable with name " << matvar->name << std::endl;

            int nn = matvar->dims[0];
            int pp = matvar->dims[1];

            init(nn,pp);

            sparse_t* sparse = static_cast<sparse_t*>(matvar->data);
            _nz = sparse->nzmax;

            double *data = static_cast<double*>(sparse->data);

            idxType current_col = 0;
            for(size_t k = 0; k < _nz; ++k)
            {
                idxType i = sparse->ir[k];
                valType val = data[k];
                if(k < sparse->jc[sparse->njc-1])
                {
                    assert(current_col < sparse->njc);
                    while(sparse->jc[current_col + 1] <= k) // look for the last idx of jc such that jc[idx+1] > k
                    {
                        current_col++;
                    }
                    idxType j = current_col;
                    (*this)(i,j)=val;
                }
            }
            matvar->mem_conserve = 1;
            Mat_VarFree(matvar);
            Mat_Close(mat);

            refreshNZ();
        }
    #else
        std::cerr << "You have to compile OpenMEEG with MATIO to read matlab files" << std::endl;
    #endif
    }

    inline void load( const char *filename )
    {
        char extension[128];
        getNameExtension(filename,extension);
        if(!strcmp(extension,"bin") || !strcmp(extension,"BIN")) load(filename,'b');
        else if(!strcmp(extension,"txt") || !strcmp(extension,"TXT")) load(filename,'t');
        else if(!strcmp(extension,"mat") || !strcmp(extension,"MAT")) loadMat(filename);
        else {
            std::cout << "Warning : Unknown file extension " << std::endl;
            std::cout << "Assuming ASCII format " << std::endl;
            load(filename,'t');
        }
    }

    inline void load( const char *filename, const char format)
    {
        if(_RowEntry!=0) destroy();

        idxType nn,pp,ii,jj;
        valType val;

        FILE *file;
        file = fopen(filename,"r");
        if (file==NULL) {std::cerr << "Error opening file: " << filename << std::endl; exit(1);}

        // read matrix size
        readDim(file,nn,pp,format);

        cellType **RowPrev = new cellType*[nn]; for(idxType i=0;i<nn;i++) RowPrev[i]=0;
        cellType *ColPrev = 0;
        cellType *current = 0;

        init(nn,pp);

        _nz = 0;
        int nread = readCell(file,ii,jj,val,format);

        assert(nread > 0); // matrix is not empty

        while (nread > 0)
        {
            assert(ii < _n);
            assert(jj < _p);

            _nz++;
            current = new cellType;
            current->i = ii;
            current->j = jj;
            current->val = val;

            if (ColPrev == 0) {
                _ColEntry[jj] = current;
                current->up = 0;
            } else {
                ColPrev->down = current;
                current->up = ColPrev;
            }

            ColPrev = current;

            if (RowPrev[ii] == 0) {
                _RowEntry[ii] = current;
                current->left = 0;
            } else {
                RowPrev[ii]->right = current;
                current->left = RowPrev[ii];
            }

            RowPrev[ii] = current;

            idxType oldjj = jj;
            nread = readCell(file,ii,jj,val,format);

            // if start new column
            if(jj!=oldjj) {current->down=0; ColPrev=0;}
        }

        // Update last cell of the matrix
        current->down = 0;
        for (idxType i=0;i<_n;i++)
            if (RowPrev[i]!=0) RowPrev[i]->right=0;

        fclose(file);

        refreshNZ();

        delete[] RowPrev;
    }

    inline void save( const char *filename ) const
    {
        char extension[128];
        getNameExtension(filename,extension);
        if(!strcmp(extension,"bin") || !strcmp(extension,"BIN")) saveBin(filename);
        else if(!strcmp(extension,"txt") || !strcmp(extension,"TXT")) saveTxt(filename);
        else if(!strcmp(extension,"mat") || !strcmp(extension,"MAT")) saveMat(filename);
        else {
            std::cout << "Warning : Unknown file extension " << std::endl;
            std::cout << "Assuming ASCII format " << std::endl;
            saveTxt(filename);
        }
    }

    inline void write(std::ostream& f) const
    {
        ((TsparseMatrix<T,I>*)this)->refreshNZ();
        f.write(reinterpret_cast<const char*>(&_n),(std::streamsize)sizeof(idxType)); // nb lines
        f.write(reinterpret_cast<const char*>(&_p),(std::streamsize)sizeof(idxType)); // nb cols

        TsparseMatrixIterator<T,I> it(*this);
        for(it.begin(); !it.end(); it.next()) {
            f.write(reinterpret_cast<const char*>(&(it.current()->i)),(std::streamsize)sizeof(idxType));
            f.write(reinterpret_cast<const char*>(&(it.current()->j)),(std::streamsize)sizeof(idxType));
            f.write(reinterpret_cast<const char*>(&(it.current()->val)),(std::streamsize)sizeof(valType));
        }
    }

    inline void read(std::istream& f)
    {
        if(_RowEntry!=0) destroy();
        f.seekg (0, std::ios::beg);
        f.read(reinterpret_cast<char*>(&_n),(std::streamsize)sizeof(idxType)); // nb lines
        f.read(reinterpret_cast<char*>(&_p),(std::streamsize)sizeof(idxType)); // nb cols
        init(_n,_p);
        idxType i,j;
        valType val;
        while (! f.eof() ) {
            f.read(reinterpret_cast<char*>(&i),(std::streamsize)sizeof(idxType));
            f.read(reinterpret_cast<char*>(&j),(std::streamsize)sizeof(idxType));
            f.read(reinterpret_cast<char*>(&val),(std::streamsize)sizeof(valType));
            (*this)(i,j)=val;
        }
        refreshNZ();
    }

    inline TsparseMatrix<T,I> operator+( TsparseMatrix<T,I>& opd)
    {
        refreshNZ();
        opd.refreshNZ();

        if(_p!=opd.ncol() || _n!=opd.nlin()) {std::cerr<<"Erreur Somme Matrice"; exit(1);}

        TsparseMatrix<T,I> result(_n,_p);

        for(idxType i=0;i<_n;i++)
        {
            cellType* startRowD=(opd._RowEntry)[i];
            cellType* startRowG=_RowEntry[i];
            if(startRowD!=0 || startRowG!=0)
            {
                while(startRowD!=0 || startRowG!=0)
                {
                    if(startRowD!=0 && startRowG!=0)
                    {
                        if(startRowD->j > startRowG->j) {(result)(i,startRowG->j)=startRowG->val; startRowG=startRowG->right;}
                        else if(startRowD->j == startRowG->j) {(result)(i,startRowD->j)=startRowD->val+startRowG->val; startRowD=startRowD->right; startRowG=startRowG->right;}
                        else if(startRowD->j < startRowG->j) {(result)(i,startRowD->j)=startRowD->val; startRowD=startRowD->right;}
                    }
                    if(startRowD==0 && startRowG!=0) {(result)(i,startRowG->j)=startRowG->val; startRowG=startRowG->right;}
                    if(startRowG==0 && startRowD!=0) {(result)(i,startRowD->j)=startRowD->val; startRowD=startRowD->right;}
                }
            }
        }

        result.refreshNZ();
        return result;
    }

    inline TsparseMatrix<T,I>& operator+=( TsparseMatrix<T,I>& opd)
    {
        refreshNZ();
        opd.refreshNZ();

        if(_p!=opd.ncol() || _n!=opd.nlin()) {std::cerr<<"Bad dimensions for sparse matrice addition"; exit(1);}

        for(idxType i=0;i<_n;i++)
        {
            cellType* startRowD=(opd._RowEntry)[i];
            cellType* startRowG=_RowEntry[i];
            if(startRowD!=0 || startRowG!=0)
            {
                while(startRowD!=0 || startRowG!=0)
                {
                    if(startRowD!=0 && startRowG!=0)
                    {
                        if(startRowD->j > startRowG->j) {(*this)(i,startRowG->j)=startRowG->val; startRowG=startRowG->right;}
                        else if(startRowD->j == startRowG->j) {(*this)(i,startRowD->j)=startRowD->val+startRowG->val; startRowD=startRowD->right; startRowG=startRowG->right;}
                        else if(startRowD->j < startRowG->j) {(*this)(i,startRowD->j)=startRowD->val; startRowD=startRowD->right;}
                    }
                    if(startRowD==0 && startRowG!=0) {(*this)(i,startRowG->j)=startRowG->val; startRowG=startRowG->right;}
                    if(startRowG==0 && startRowD!=0) {(*this)(i,startRowD->j)=startRowD->val; startRowD=startRowD->right;}
                }
            }
        }
        return *this;
    }

    inline TsparseMatrix<T,I> operator-()
    {
        TsparseMatrix<T,I> result=(*this);

        for(idxType i=0;i<_n;i++)
        {
            cellType *pt;//,*ptbis;
            pt=result._RowEntry[i];
            while(pt!=0)
            {
                pt->val=-pt->val;
                pt=pt->right;
            }
        }
        result.refreshNZ();
        return result;
    }

    inline TsparseMatrix<T,I> operator*(T alpha)
    {
        TsparseMatrix<T,I> result=(*this);

        for(idxType i=0;i<_n;i++)
        {
            cellType *pt;//,*ptbis;
            pt=result._RowEntry[i];
            while(pt!=0)
            {
                pt->val=alpha*pt->val;
                pt=pt->right;
            }
        }
        result.refreshNZ();
        return result;
    }

    inline TsparseMatrix<T,I> transpose() const
    {
        TsparseMatrix<T,I> result=(*this);
        result.autotranspose();
        return result;
    }

    inline void autotranspose()
    {
        cellType **tempCell;
        tempCell=_RowEntry;
        _RowEntry=_ColEntry;
        _ColEntry=tempCell;
        idxType  tempInt;
        tempInt=_n;
        _n=_p;
        _p=tempInt;

        tempInt=_nnzrows;
        _nnzrows=_nnzcols;
        _nnzcols=tempInt;

        idxType *tempIntp;
        tempIntp=_nzrows;
        _nzrows=_nzcols;
        _nzcols=tempIntp;

        if(_sp==row) _sp=column;
        else _sp=row;

        for(idxType i=0;i<_p;i++)
        {
            cellType *pt;//*ptbis;
            pt=_ColEntry[i];
            while(pt!=0)
            {
                cellType *tempPt;
                tempPt=pt->left;
                pt->left=pt->up;
                pt->up=tempPt;
                tempPt=pt->right;
                pt->right=pt->down;
                pt->down=tempPt;
                tempInt=pt->j;
                pt->j=i;
                pt->i=tempInt;
                pt=pt->down;
            }
        }
    }

    inline valType operator()(size_t i, size_t j) const
    {
        cellType **Entry,*temp,ruse;
        valType result=0;

        assert(i<_n && j<_p);

        if (_sp==row) Entry=_RowEntry; else Entry=_ColEntry;

        if (_RowEntry[i]==0 || _ColEntry[j]==0) result=0;
        else {
            bool status=false;

            if (_sp==row)
            {
                ruse.right=_RowEntry[i];
                temp=&ruse;
                do
                {
                    temp=temp->right;
                    if(temp->j==j) {result=temp->val; status=true;} //On a trouve, on s'arrete
                    if(temp->j>j || temp->right==0) status=true;    //On n'a pas trouve, on s'arrete
                }
                while(!status);
            }
            else
            {
                ruse.down=_ColEntry[j];
                temp=&ruse;
                do
                {
                    temp=temp->down;
                    if(temp->i==i) {result=temp->val; status=true;} //On a trouve, on s'arrete
                    if(temp->i>i || temp->down==0) status=true;     //On n'a pas trouve, on s'arrete
                }
                while(!status);
            }
        }
        return result;

    }; //end operator ()

    inline T& operator()(size_t i, size_t j)
    {
        assert(i<_n && j<_p);

        cellType **Entry,*temprow=0,*tempcol=0,*temp=0,ruse,*insRow=0,*insCol=0;
        bool flag=false,status=false;
        ruse.up=0;
        ruse.down=0;
        ruse.left=0;
        ruse.right=0;
        temp=&ruse;

        if(_sp==row) Entry=_RowEntry; else Entry=_ColEntry;

        if(_RowEntry[i]!=0 && _ColEntry[j]!=0)
        {
            status=false;

            if(_sp==row)
            {
                ruse.right=_RowEntry[i];
                temprow=&ruse;
                do
                {
                    temprow=temprow->right;
                    if(temprow->j==j) {return(temprow->val);}       //On a trouve, on s'arrete
                    if(temprow->j>j) {flag=true; status=true;}   //On a trouve un elt d'index >, on s'arrete
                    if(temprow->right==0) status=true;              //On n'a pas trouve avant la fin
                }
                while(!status);
            }
            else
            {
                ruse.down=_ColEntry[j];
                tempcol=&ruse;
                do
                {
                    tempcol=tempcol->down;
                    if(tempcol->i==i) {return(tempcol->val); } //On a trouve, on s'arrete
                    if(tempcol->i>i) {flag=true; status=true;} //On a trouve un elt d'index >, on s'arrete
                    if(tempcol->down==0) status=true;   //On n'a pas trouveavant la fin
                }
                while(!status);
            }

        }

        // si on arrive à ce stade, c'est que l'element cherche n'existe pas, il faut donc le creer

        // On cherche l'endroit ou l'insÈrer suivant la ligne concernÈe
        if(_RowEntry[i]!=0)
        {
            status=false;
            flag=false;
            ruse.right=_RowEntry[i];
            temprow=&ruse;
            do
            {
                temprow=temprow->right;
                if(temprow->j>j) {flag=true; status=true;}   //On a trouve un elt d'index >, on s'arrete
                if(temprow->right==0) status=true;              //On n'a pas trouve avant la fin
            }
            while(!status);

            if(flag) insRow=temprow->left;  //element apres lequel doit Ítre insÈrÈ le nouveau
            else insRow=temprow;            //on indique qu'il faut l'ajouter à la fin
        }

        // On cherche l'endroit ou l'insÈrer suivant la colonne concernÈe
        if(_ColEntry[j]!=0)
        {
            status=false;
            flag=false;
            ruse.down=_ColEntry[j];
            tempcol=&ruse;
            do
            {
                tempcol=tempcol->down;
                if(tempcol->i>i) {flag=true; status=true;}  //On a trouve un elt d'index >, on s'arrete
                if(tempcol->down==0) status=true;           //On n'a pas trouveavant la fin
            }
            while(!status);

            if(flag) insCol=tempcol->up;    //element apres lequel doit Ítre insÈrÈ le nouveau
            else insCol=tempcol;            //on indique qu'il faut l'ajouter à la fin
        }


        cellType *tempbis=new cellType(); _nz++; //Un element non nul de plus dans la matrice
        tempbis->i=i;
        tempbis->j=j;
        //result=temp->val;

        if(insRow==0 || _RowEntry[i]==0) // Si l'on ajoute une cellule en tete de ligne
        {
            tempbis->left=0;

            if(_RowEntry[i]==0) // l'element ajouter est seul sur cette ligne
            {
                tempbis->right=0;
            }
            else
            {
                tempbis->right=_RowEntry[i];//insRow->right;
                _RowEntry[i]->left=tempbis;
            }

            _RowEntry[i]=tempbis;
        }
        else //on doit ajouter ailleurs qu'en tete  de ligne
        {
            tempbis->left=insRow;
            tempbis->right=insRow->right;
            insRow->right=tempbis;
            if(tempbis->right!=0) tempbis->right->left=tempbis;
        }

        if(insCol==0 || _ColEntry[j]==0) // Si l'on ajoute une cellule en tete de colonne
        {
            tempbis->up=0;

            if(_ColEntry[j]==0) // l'element ajouter est seul sur cette colonne
            {
                tempbis->down=0;
            }
            else
            {
                tempbis->down=_ColEntry[j];//insCol->down;
                _ColEntry[j]->up=tempbis;
            }

            _ColEntry[j]=tempbis;
        }
        else //on doit ajouter ailleurs qu'en tete de colonne
        {
            tempbis->up=insCol;
            tempbis->down=insCol->down;
            insCol->down=tempbis;
            if(tempbis->down!=0) tempbis->down->up=tempbis;

        }

        return tempbis->val;
    }; //fin de operator ()

    inline const TsparseMatrix<T,I>& operator=(const TsparseMatrix<T,I>& modele)
    {
        // Delete current instance
        if(_RowEntry!=0) destroy();

        // Recreate new content
        _n = modele._n;
        _p = modele._p;
        _nz = modele._nz;
        _sp = modele._sp;
        _tr = modele._tr;
        _RowEntry = new cellType*[_n]; for(idxType i=0;i<_n;i++) _RowEntry[i]=0;
        _ColEntry = new cellType*[_p]; for(idxType j=0;j<_p;j++) _ColEntry[j]=0;
        _nzrows = new idxType[_n]; for(idxType i=0;i<_n;i++) _nzrows[i]=modele._nzrows[i];
        _nzcols = new idxType[_p]; for(idxType i=0;i<_p;i++) _nzcols[i]=modele._nzcols[i];
        _nnzrows = modele._nnzrows;
        _nnzcols = modele._nnzcols;

        for(idxType i=0;i<_n;i++)
        {
            cellType *pt;
            pt=modele._RowEntry[i];
            while(pt!=0)
            {
                (*this)(pt->i,pt->j) = modele(pt->i,pt->j);
                pt=pt->right;
            }
        }
        return *this;
    }
};

template<class T,class I> class TsparseMatrixIterator
{
    public:

        typedef MyCell<T,I> cellType;

        inline TsparseMatrixIterator(const TsparseMatrix<T,I> &M) : _M(M) {
            _current = _M.first();
        }

        inline void begin(void) {
            _current = _M.first();
        }

        inline void next(void) {
            if (_current != NULL) {
                I next_col = _current->j + 1;
                _current = _current->down;
                while (_current == NULL && next_col < _M.ncol()) {
                    _current = _M.getColEntry()[next_col];
                    next_col++;
                }
            }
        }

        inline int end(void) { return (_current == NULL); }

        inline cellType *current(void) {
            if(_current == NULL) return NULL; else return _current;
        }

    private:
        cellType *_current;
        const TsparseMatrix<T,I>& _M;
};

template <> inline vecteur sparse_matrice::operator*( const vecteur &x) const
{
    vecteur ret(nlin());
    ret.set(0);

    for(idxType i=0;i<_nnzrows;i++)
    {
        MyCell<sparse_matrice::valType,sparse_matrice::idxType>* startRow=_RowEntry[_nzrows[i]];
        if(startRow!=0)
        {
            double total=0;
            while(startRow!=0)
            {
                total+=startRow->val*x(startRow->j);
                startRow=startRow->right;
            }
            ret(_nzrows[i])=total;
        }
    }

    return ret;
}

#endif
