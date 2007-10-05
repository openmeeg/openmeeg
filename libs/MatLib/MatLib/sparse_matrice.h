#ifndef H_SPARSE_MATRICE
#define H_SPARSE_MATRICE

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <fstream>

#include "vecteur.h"
#include "generic_matrix.h"

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

template<class T,class I> class TsparseMatrix: public TgenericMatrix<T>
{
public:
    typedef T valType;
    typedef I idxType;
    typedef MyCell<T,I> cellType;
protected:
    // On a 2 tableaux d'indexs pour parcourir les matrices
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
    inline cellType **getColEntry() {return _ColEntry;}
    inline cellType **getRowEntry() {return _RowEntry;}
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

    inline TsparseMatrix ( const TsparseMatrix<T,I>& modele )
    {

        _n=modele._n;
        _p=modele._p;
        //nz=modele.nz;
        _nz=0;
        _sp=modele._sp;
        _tr=modele._tr;
        _RowEntry=new cellType*[_n]; for(idxType i=0;i<_n;i++) _RowEntry[i]=0;
        _ColEntry=new cellType*[_p]; for(idxType j=0;j<_p;j++) _ColEntry[j]=0;
        _nzrows=new idxType[_n]; for(idxType i=0;i<_n;i++) _nzrows[i]=modele._nzrows[i];
        _nzcols=new idxType[_p]; for(idxType i=0;i<_p;i++) _nzcols[i]=modele._nzcols[i];
        _nnzrows=modele._nnzrows;
        _nnzcols=modele._nnzcols;

        for(idxType i=0;i<_n;i++)
        {
            cellType *pt;
            pt=modele._RowEntry[i];
            while(pt!=0)
            {
                (*this)(pt->i,pt->j)=modele(pt->i,pt->j);
                pt=pt->right;
            }
        }
    }

    vecteur operator * ( const vecteur &x) const;

    inline TsparseMatrix<T,I> operator * (  TsparseMatrix<T,I>& opd)
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

    inline void saveBin( const char * filename) const // Les colonnes sont ecrites les unes à la suite des autres
    {
        FILE *file;
        file=fopen(filename,"wb");
	    if(file==NULL) {std::cerr<<"Error opening file: " << filename << std::endl; exit(1);}

        // On ecrit avant toute chose le nombre de lignes et de colonnes de la matrice
        idxType nn=this->_n;
        idxType pp=this->_p;
        fwrite(&nn,sizeof(idxType),1,file);
        fwrite(&pp,sizeof(idxType),1,file);

        // On parcourt les colonnes
        for(idxType j=0;j<_p;j++)
        {
            if(_ColEntry[j]!=0) // Si la colonne n'est pas vide, il y a quelque chos à ecrire
            {
                cellType *ce=_ColEntry[j];
                while(ce!=0)
                {
                    fwrite(&(ce->i),sizeof(idxType),1,file);
                    fwrite(&(ce->j),sizeof(idxType),1,file);
                    fwrite(&(ce->val),sizeof(T),1,file);
                    ce=ce->down;
                }
            }
        } //fin de la boucle sur les colonnes

        fclose(file);
    }

    inline void saveTxt( const char * filename) const // Les colonnes sont ecrites les unes à la suite des autres
    {
        std::ofstream of(filename);

        of<<(unsigned long int)_n<<"\t"<<(unsigned long int)_p<<std::endl;
        of<<(*this);
    }

    inline void loadTxt( const char * filename) // Les colonnes sont ecrites les unes à la suite des autres
    {
        if(_RowEntry!=0) destroy();

        idxType nn,pp,ii,jj;
        T vval;

        FILE *file;
        file=fopen(filename,"r");
	    if(file==NULL) {std::cerr<<"Error opening file: " << filename << std::endl; exit(1);}

        // On ecrit avant toute chose le nombre de lignes et de colonnes de la matrice
        int nnn,ppp,iii,jjj;
        fscanf(file,"%d",&nnn);
        fscanf(file,"%d",&ppp);
        nn=(size_t)nnn; pp=(size_t)ppp;

        cellType **RowPrev=new cellType*[nn]; for(idxType i=0;i<nn;i++) RowPrev[i]=0;
        cellType *ColPrev=0;
        cellType *current=0;

        init(nn,pp);

        // On parcourt le fichier et on remplit la matrice
        _nz=0;
        int nread=1;
        fscanf(file,"%d",&iii);
        fscanf(file,"%d",&jjj);
        fscanf(file,"%lf",&vval);
        ii=iii;jj=jjj;
        while(nread>0)
        {
            _nz++;

            current=new cellType;
            current->i=ii;
            current->j=jj;
            current->val=vval;

            // On met le pointeur up àjour
            if(ColPrev==0)
            {
                _ColEntry[jj]=current;
                current->up=0;
            }
            else
            {
                ColPrev->down=current;
                current->up=ColPrev;
            }

            ColPrev=current; // On met à jour le ColPrev

            //On met à jour le pointeur left
            if(RowPrev[ii]==0)
            {
                _RowEntry[ii]=current;
                current->left=0;
            }
            else
            {
                RowPrev[ii]->right=current;
                current->left=RowPrev[ii];
            }

            RowPrev[ii]=current; // On met à jour le RowPrev[ii]
            ColPrev=current;

            // On lit les valeurs suivantes dans le fichier
            nread=0;
            idxType oldjj=jj;
            nread+=fscanf(file,"%d",&iii);
            nread+=fscanf(file,"%d",&jjj);
            ii=iii; jj=jjj;
            nread+=fscanf(file,"%lf",&vval);

            // Si on commence une nouvelle colonne, on ferme l'ancienne
            if(jj!=oldjj) {current->down=0; ColPrev=0;}
        }

        // Une fois qu'on à fini la lecture de la matrice, il convient 
        // de verouiller la fin de colonne et les fins de lignes
        current->down=0;
        for(idxType i=0;i<_n;i++)
            if(RowPrev[i]!=0) RowPrev[i]->right=0;

        // On referme le fichier
        fclose(file);

        // On met à jour les nnzs
        refreshNZ();

        //On efface les variables temporaires
        delete[] RowPrev;
    }

    inline void loadBin( const char * filename) // Les colonnes sont ecrites les unes à la suite des autres
    {
        if(_RowEntry!=0) destroy();

        idxType nn,pp,ii,jj;
        T vval;

        FILE *file;
        file=fopen(filename,"rb");
	    if(file==NULL) {std::cerr<<"Error opening sparse matrix binary file: " << filename << std::endl; exit(1);}

        // On ecrit avant toute chose le nombre de lignes et de colonnes de la matrice
        fread(&nn,sizeof(idxType),1,file);
        fread(&pp,sizeof(idxType),1,file);

        cellType **RowPrev=new cellType*[nn]; for(idxType i=0;i<nn;i++) RowPrev[i]=0;
        cellType *ColPrev=0;
        cellType *current=0;

        init(nn,pp);

        // On parcourt le fichier et on remplit la matrice
        _nz=0;
        size_t nread=1;
        fread(&ii,sizeof(idxType),1,file);
        fread(&jj,sizeof(idxType),1,file);
        fread(&vval,sizeof(T),1,file);
        while(nread>0)
        {
            _nz++;

            current=new cellType;
            current->i=ii;
            current->j=jj;
            current->val=vval;

            // On met le pointeur up à jour
            if(ColPrev==0)
            {
                _ColEntry[jj]=current;
                current->up=0;
            }
            else
            {
                ColPrev->down=current;
                current->up=ColPrev;
            }

            ColPrev=current; // On met à jour le ColPrev

            //On met à jour le pointeur left
            if(RowPrev[ii]==0)
            {
                _RowEntry[ii]=current;
                current->left=0;
            }
            else
            {
                RowPrev[ii]->right=current;
                current->left=RowPrev[ii];
            }

            RowPrev[ii]=current; // On met à jour le RowPrev[ii]
            ColPrev=current;

            // On lit les valeurs suivantes dans le fichier
            nread=0;
            idxType oldjj=jj;
            nread+=fread(&ii,sizeof(idxType),1,file);
            nread+=fread(&jj,sizeof(idxType),1,file);
            nread+=fread(&vval,sizeof(T),1,file);

            // Si on commence une nouvelle colonne, on ferme l'ancienne
            if(jj!=oldjj) {current->down=0; ColPrev=0;}
        }

        // Une fois qu'on à fini la lecture de la matrice, il convient de verouiller la fin de colonne et les fins de lignes
        current->down=0;
        for(idxType i=0;i<_n;i++)
            if(RowPrev[i]!=0) RowPrev[i]->right=0;

        // On referme le fichier
        fclose(file);

        // On met à jour les nnzs
        refreshNZ();

        //On efface les variables temporaires
        delete[] RowPrev;
    }

    inline void write(std::ostream& f) const
    {
        ((TsparseMatrix<T,I>*)this)->refreshNZ();
        f.write((const char*)&_n,(std::streamsize)sizeof(idxType)); // n lignes
        f.write((const char*)&_p,(std::streamsize)sizeof(idxType)); // n col
        f.write((const char*)&_nz,(std::streamsize)sizeof(idxType)); // n non nuls
        // On parcourt les colonnes
        for(idxType j=0;j<_p;j++)
        {
            if(_ColEntry[j]!=0) // Si la colonne n'est pas vide, il y a quelque chos à ecrire
            {
                cellType *ce=_ColEntry[j];
                while(ce!=0)
                {
                    f.write((const char*)&(ce->i),(std::streamsize)sizeof(idxType));
                    f.write((const char*)&(ce->j),(std::streamsize)sizeof(idxType));
                    f.write((const char*)&(ce->val),(std::streamsize)sizeof(T));
                    ce=ce->down;
                }
            }
        } //fin de la boucle sur les colonnes
    }

    inline void read(std::istream& f)
    {
        // On parcourt les colonnes
        destroy();
        f.read((char*)&_n,(std::streamsize)sizeof(idxType)); // n lignes
        f.read((char*)&_p,(std::streamsize)sizeof(idxType)); // n col
        init(_n,_p);
        idxType nz; // FIXME:
        f.read((char*)&nz,(std::streamsize)sizeof(idxType)); // n non nuls
        for(idxType k=0;k<nz;k++)
        {
            idxType i,j;
            T val;

            f.read((char*)&i,(std::streamsize)sizeof(idxType));
            f.read((char*)&j,(std::streamsize)sizeof(idxType));
            f.read((char*)&val,(std::streamsize)sizeof(T));
            (*this)(i,j)=val;
        }
        refreshNZ();
    }

    inline void display() const //Procedure lente car pas optimisee!!
    {
        bool flag=true;
        for(idxType i=0;i<this->_n;i++)
            for(idxType j=0;j<this->_p;j++)
        {
            double val=((const TsparseMatrix<T,I> &)(*this))(i,j);
            if(val!=0)
            {
                flag=false;
                std::cout<<"\t("<<(unsigned int)i<<","<<(unsigned int)j<<") = "<<val<<std::endl;
            }
        }

        if(flag) std::cout<<"Matrice Vide"<<std::endl;
    }

    inline TsparseMatrix<T,I> operator + ( TsparseMatrix<T,I>& opd)
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

    inline void operator += ( TsparseMatrix<T,I>& opd)
    {
        refreshNZ();
        opd.refreshNZ();

        if(_p!=opd.ncol() || _n!=opd.nlin()) {std::cerr<<"Erreur Somme Matrice"; exit(1);}

        //TsparseMatrix<T,I> *this(this->n,this->p);

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
    }

    inline TsparseMatrix<T,I> operator - ()
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

    inline TsparseMatrix<T,I> operator * (T alpha)
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

    inline TsparseMatrix<T,I> transpose()
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

    inline T operator () (size_t i, size_t j) const
    {
        cellType **Entry,*temp,ruse;
        T result=0;

        assert(i<_n && j<_p);

        if(_sp==row) Entry=_RowEntry; else Entry=_ColEntry;

        if(_RowEntry[i]==0 || _ColEntry[j]==0) result=0;
        else
        {
            bool status=false;

            if(_sp==row)
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

    }; //fin de operator ()

    inline T& operator () (size_t i, size_t j)
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

    inline void operator = (TsparseMatrix<T,I>& modele)
    {
        // On efface le contenu precedent de la matrice

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

        if (_RowEntry)
        {
            delete[] _RowEntry;
            delete[] _ColEntry;
            delete[] _nzrows;
            delete[] _nzcols;
        }

        // On recree le nouveau contenu

        _n=modele._n;
        _p=modele._p;
        _nz=modele._nz;
        _sp=modele._sp;
        _tr=modele._tr;
        _RowEntry=new cellType*[_n]; for(idxType i=0;i<_n;i++) _RowEntry[i]=0;
        _ColEntry=new cellType*[_p]; for(idxType j=0;j<_p;j++) _ColEntry[j]=0;
        _nzrows=new idxType[_n]; for(idxType i=0;i<_n;i++) _nzrows[i]=modele._nzrows[i];
        _nzcols=new idxType[_p]; for(idxType i=0;i<_p;i++) _nzcols[i]=modele._nzcols[i];
        _nnzrows=modele._nnzrows;
        _nnzcols=modele._nnzcols;

        for(idxType i=0;i<_n;i++)
        {
            cellType *pt;
            pt=modele._RowEntry[i];
            while(pt!=0)
            {
                (*this)(pt->i,pt->j)=modele(pt->i,pt->j);
                pt=pt->right;
            }
        }
    }

};

template<class T,class I> inline std::ostream& operator<<(std::ostream& f,const TsparseMatrix<T,I> &M)
{
    // Loop on columns
    for(size_t j=0;j<M.ncol();j++)
    {
        if(((TsparseMatrix<T,I> &)M).getColEntry()[j]!=0) // Si la colonne n'est pas vide, il y a quelque chos à ecrire
        {
            typename TsparseMatrix<T,I>::cellType *ce=((TsparseMatrix<T,I> &)M).getColEntry()[j];
            while(ce!=0)
            {
                f<<(unsigned long int)ce->i<<"\t"<<(unsigned long int)ce->j<<"\t"<<ce->val<<std::endl;
                ce=ce->down;
            }
        }
    }
    return f;
}

typedef TsparseMatrix<double,size_t> sparse_matrice;

template <> inline vecteur sparse_matrice::operator * ( const vecteur &x) const
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

template class TsparseMatrix<double,size_t>;

// template inline std::ostream& operator<<(std::ostream& f,const TsparseMatrix<double,size_t> &M);

#endif
