#ifndef H_sparse_matrice
#define H_sparse_matrice

#include "vecteur.h"
#include "MySparseLib.h"

typedef MySparseLib<double,size_t> sparse_matrice;

template class MySparseLib<double,size_t>;

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

template inline std::ostream& operator<<(std::ostream& f,const MySparseLib<double,size_t> &M);


#endif
