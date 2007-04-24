#ifndef H_sparse_vecteur
#define H_sparse_vecteur

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <fstream>
#include <list>
#include "vecteur.h"

template<class T,class I> struct VCell
{
    T val;
    I i;
    VCell<T,I> *next;
    VCell() {next=0;}
};


template<class T,class I> class sparse_vecteur
{

public:
    typedef T valType;
    typedef I idxType;
    typedef VCell<T,I> cellType;
    typedef sparse_vecteur<T,I> self;

    class iterateur {
        cellType *pos;
    public:
        iterateur(cellType *pos) {this->pos=pos;}
        inline idxType i() const {return pos->i;}
        inline valType val() const {return pos->val;}
        inline const iterateur& operator++() {pos=pos->next; return *this;}
        inline bool operator== (const iterateur& it) {return(pos==it.pos);}
        inline bool operator!= (const iterateur& it) {return(pos!=it.pos);}
    };
    iterateur begin() const {return iterateur(Entry);}
    iterateur end() const {return iterateur(0);}

protected:
    cellType *Entry;
    idxType n;
    idxType nz;
    int *count;

    void alloc(size_t N) {n=N; nz=0; count=new int; *count=1; Entry=0;}
    void copy(const self& b);
    void destroy();
    // insere une nouvelle cellule avant / detruit **cptr (cptr doit etre &Entry ou un &(VCell::next))
    void insere(cellType **cptr);
    void erase(cellType **cptr);

public:
    sparse_vecteur(size_t N=0) {alloc(N);}
    sparse_vecteur(const sparse_vecteur& v) {copy(v);}
    ~sparse_vecteur() {destroy();}
    const self& operator= (const self& v) {destroy(); copy(v); return v;}

    inline size_t size() const {return (size_t)n;}
    inline size_t getNz() const {return (size_t)nz;}
    self duplicate() const;
    void operator+= (const self& b);
    void operator-= (const self& b);
    void operator*= (T x);
    void operator/= (T x);
    self operator+ (const self& b) const {self c=duplicate(); c+=b; return c;} // pas optimal
    self operator- (const self& b) const {self c=duplicate(); c-=b; return c;} // pas optimal
    self operator* (T x) const {self c=duplicate(); c*=x; return c;} // pas optimal
    self operator/ (T x) const {self c=duplicate(); c/=x; return c;} // pas optimal
    void addmult(double k, const self& b, double tol=0); // realise : *this += k*b

    T operator() (idxType i) const;
    T& operator() (idxType i);
    cellType* getEntry() {return Entry;}

    friend std::ostream& operator<< (std::ostream& out, const self& v) {
        for (cellType* c=v.Entry; c!=0; c=c->next) std::cout<<c->i<<":"<<c->val<<" ";
        return out;
    }

};

template<class T,class I> void sparse_vecteur<T,I>::copy(const sparse_vecteur<T,I>& v)
{
    count=v.count;
    (*count)++;
    n=v.n;
    nz=v.nz;
    Entry=v.Entry;
}

template<class T,class I> void sparse_vecteur<T,I>::destroy()
{
    (*count)--;
    if (*count==0) {
        delete count;
        cellType *pt,*ptbis;
        pt=Entry;
        while(pt!=0)
        {
            ptbis=pt->next;
            delete pt;
            pt=ptbis;
        }
    }
}

template<class T,class I> void sparse_vecteur<T,I>::insere(VCell<T,I> **cptr)
{
    cellType* newcell=new cellType;
    newcell->next=(*cptr);
    *cptr=newcell;
    nz++;
}

template<class T,class I> void sparse_vecteur<T,I>::erase(VCell<T,I> **cptr)
{
    cellType* oldcell=*cptr;
    *cptr=oldcell->next;
    delete oldcell;
    nz--;
}

template<class T,class I> sparse_vecteur<T,I> sparse_vecteur<T,I>::duplicate() const
{
    self b;
    b.alloc(n);
    cellType *c=Entry, **bcptr=&(b.Entry);
    while(c!=0) {
        *bcptr = new cellType;
        (*bcptr)->i=c->i;
        (*bcptr)->val=c->val;
        c = c->next;
        bcptr = &((*bcptr)->next);
    }
    return b;
}

template<class T,class I> void sparse_vecteur<T,I>::operator+= (const self& b)
{
    n = max(n,b.n);
    cellType **cptr=&Entry, *bc=b.Entry;
    while(bc!=0) {
        while((*cptr)!=0 && (*cptr)->i<bc->i) {cptr=&((*cptr)->next);}
        if ((*cptr)!=0 && (*cptr)->i==bc->i)
        {
            (*cptr)->val+=bc->val;
            if ((*cptr)->val==0) erase(cptr);
        }
        else
        {
            insere(cptr);
            (*cptr)->i=bc->i;
            (*cptr)->val=bc->val;
        }
        bc=bc->next;
    }
}

template<class T,class I> void sparse_vecteur<T,I>::operator-= (const self& b)
{
    n = max(n,b.n);
    cellType **cptr=&Entry, *bc=b.Entry;
    while(bc!=0)
    {
        while((*cptr)!=0 && (*cptr)->i<bc->i) {cptr=&((*cptr)->next);}
        if ((*cptr)!=0 && (*cptr)->i==bc->i)
        {
            (*cptr)->val-=bc->val;
            if ((*cptr)->val==0) erase(cptr);
        }
        else
        {
            insere(cptr);
            (*cptr)->i=bc->i;
            (*cptr)->val=-bc->val;
        }
        bc=bc->next;
    }
}

template<class T,class I> void sparse_vecteur<T,I>::operator*= (T x)
{
    if (x==0) {destroy(); alloc(0);}
    for (cellType *c=Entry; c!=0; c=c->next) c->val*=x;
}

template<class T,class I> void sparse_vecteur<T,I>::operator/= (T x)
{
    assert(x!=0);
    for (cellType *c=Entry; c!=0; c=c->next) c->val/=x;
}

template<class T,class I> void sparse_vecteur<T,I>::addmult(double k, const self& b, double tol)
{
    n = max(n,b.n);
    cellType **cptr=&Entry, *bc=b.Entry;
    while(bc!=0)
    {
        while((*cptr)!=0 && (*cptr)->i<bc->i) {cptr=&((*cptr)->next);}
        if ((*cptr)!=0 && (*cptr)->i==bc->i)
        {
            (*cptr)->val+=k*bc->val;
            if (abs((*cptr)->val)<=tol) erase(cptr);
        }
        else if (abs(k*bc->val)>tol)
        {
            insere(cptr);
            (*cptr)->i=bc->i;
            (*cptr)->val=k*bc->val;
        }
        bc=bc->next;
    }
}

template<class T,class I> T sparse_vecteur<T,I>::operator() (I i) const
{
    cellType* c;
    for (c=Entry; c!=0 && c->i<i; c=v->next) {};
    return (c!=0 && c->i==i)?c->val:0;
}

template<class T,class I> T& sparse_vecteur<T,I>::operator() (I i)
{
    cellType** cptr=&Entry; while ((*cptr)!=0 && (*cptr)->i<i) {cptr=&(*cptr)->next;}
    if ((*cptr)!=0 && (*cptr)->i==i) return (*cptr)->val;
    // si on atteint ce code, (*cptr)'est qu'il faut créer une nouvelle cellule avant *(*cptr)
    insere(cptr);
    (*cptr)->i=i;
    return (*cptr)->val;
}

//template<class T,class I> std::ostream& operator<< (std::ostream& out, const sparse_vecteur<T,I>& v)

#endif
