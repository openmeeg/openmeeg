/* FILE: $Id$ */

/*
Project Name : OpenMEEG

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

#include "MatLibConfig.h"
#include "vector.h"
#include "matrix.h"
#include "symmatrix.h"

Vector::Vector() : LinOp(0,1,FULL,ONE),t(0),count(0) {}
Vector::Vector(size_t N) : LinOp(N,1,FULL,ONE),t(0),count(0) { alloc_data(); }
Vector::Vector(double* T, int* COUNT, size_t N) : LinOp(N,1,FULL,ONE),t(T),count(COUNT) {(*count)++;}
size_t Vector::size() const { return nlin(); }
bool Vector::empty() const {return t==0;}
double*  Vector::data() const {return t;}
int* Vector::DangerousGetCount () const {return count;}
Vector Vector::duplicate() const
{
    Vector A;
    if (t)
    {
        A.nlin() = nlin();
        A.alloc_data();
        copyout(A.t);
    }
    return A;
}
void Vector::copyin(const Vector& A)
{
    assert(nlin()==A.nlin());
    copyin(A.t);
}

void Vector::operator/=(double x) {(*this)*=(1.0/x);}
Vector Vector::operator/(double x) const {return (*this)*(1.0/x);}
double Vector::mean() const {return sum()/size();}
Vector operator * (const double &d, const Vector &v) {return v*d;}

std::ostream& operator<<(std::ostream& f,const Vector &M) {
    for (size_t i=0;i<M.size();i++)
    {
        f << M(i) << " ";
    }
    return f;
}

std::istream& operator>>(std::istream& f,Vector &M) {
    for (size_t i=0;i<M.size();i++)
    {
        f >> M(i);
    }

    return f;
}

void Vector::alloc_data() {
    t=new double[nlin()];
    count=new int[1];
    (*count)=1;
}

void Vector::destroy() {
    if (t!=0)
    {
        (*count)--;
        if (*count == 0)
        {
            delete[] t;
            delete[] count;
        }
    }
}

void Vector::copy(const Vector& A) {
    t=A.t;
    nlin()=A.nlin();
    if (t)
    {
        count=A.count;
        (*count)++;
    }
}

void Vector::copyout(double * p) const {
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)nlin(),t,1,p,1);
#else
    for( size_t i=0; i<nlin(); i++ )
    {
        p[i]=t[i];
    }
#endif
}

void Vector::copyin(const double * p) {
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)nlin(),p,1,t,1);
#else
    for( size_t i=0; i<nlin(); i++ )
    {
        t[i]=p[i];
    }
#endif
}

const Vector& Vector::operator=(const Vector& A)
{
    destroy();
    copy(A);
    return *this;
}

Vector::Vector(const Vector& A)
{
    copy(A);
}

Vector::Vector(Matrix& A)
{
    nlin()=A.nlin()*A.ncol();
    t=A.data();
    if (t)
    {
        count=A.DangerousGetCount();
        (*count)++;
    }
}

Vector::Vector(SymMatrix& A)
{
    nlin()=A.nlin()*(A.nlin()+1)/2;
    t=A.data();
    if (t)
    {
        count=A.DangerousGetCount();
        (*count)++;
    }
}

Vector Vector::operator+(const Vector& v) const {
    assert(nlin()==v.nlin());
    Vector p=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)nlin(),1,v.t,1,p.t,1);
#else
    for( size_t i=0; i<nlin(); i++ )
    {
        p.t[i]=t[i]+v.t[i];
    }
#endif
    return p;
}

Vector Vector::operator-(const Vector& v) const {
    assert(nlin()==v.nlin());
    Vector p=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)nlin(),-1,v.t,1,p.t,1);
#else
    for( size_t i=0; i<nlin(); i++ )
    {
        p.t[i]=t[i]-v.t[i];
    }
#endif
    return p;
}

void Vector::operator+=(const Vector& v) {
    assert(nlin()==v.nlin());
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)nlin(),1,v.t,1,t,1);
#else
    for( size_t i=0; i<nlin(); i++ )
    {
        t[i]+=v.t[i];
    }

#endif
}

void Vector::operator-=(const Vector& v) {
    assert(nlin()==v.nlin());
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)nlin(),-1,v.t,1,t,1);
#else
    for( size_t i=0; i<nlin(); i++ )
    {
        t[i]-=v.t[i];
    }
#endif
}

double Vector::operator*(const Vector& v) const {
    assert(nlin()==v.nlin());
#ifdef HAVE_BLAS
    return BLAS(ddot,DDOT)((int)nlin(),t,1,v.t,1);
#else
    double s=0;
    for( size_t i=0; i<nlin(); i++ )
    {
        s+=t[i]*v.t[i];
    }
    return s;
#endif
}

Vector Vector::kmult(const Vector& v) const { // Kronecker multiplication
    assert(nlin() == v.nlin());
// #ifdef HAVE_BLAS
//     // FIXME : add blas version
// #else
    Vector p(nlin());
    for( size_t i=0; i<nlin(); i++ )
    {
        p(i) = v(i)*t[i];
    }
// #endif
    return p;
}

Vector Vector::operator*(double x) const {
#ifdef HAVE_BLAS
    Vector p=duplicate();
    BLAS(dscal,DSCAL)((int)nlin(),x,p.t,1);
#else
    Vector p(nlin());
    for( size_t i=0; i<nlin(); i++ )
    {
        p.t[i]=x*t[i];
    }
#endif
    return p;
}

Vector Vector::operator+(double x) const
{
    Vector p=duplicate();
    for( size_t i=0; i<nlin(); i++ )
    {
        p.t[i]+=x;
    }
    return p;
}

Vector Vector::operator-(double x) const
{
    Vector p=duplicate();
    for( size_t i=0; i<nlin(); i++ )
        p.t[i]-=x;

    return p;
}

void Vector::operator*=(double x) {
#ifdef HAVE_BLAS
    BLAS(dscal,DSCAL)((int)nlin(),x,t,1);
#else
    for( size_t i=0; i<nlin(); i++ )
    {
        t[i]*=x;
    }
#endif
}

double Vector::norm() const
{
#ifdef HAVE_BLAS
    return BLAS(dnrm2,DNRM2)((int)nlin(),t,1);
#else
    std::cout << "'Vector::norm' not implemented" << std::endl;
    exit(1);
    return 0;
#endif
}

void Vector::set(double x) {
    assert(nlin()>0);
    for( size_t i=0; i<nlin(); i++ )
    {
        t[i]=x;
    }
}

double Vector::sum() const
{
    double s=0;
    for (size_t i=0; i<nlin(); i++)
    {
        s+=t[i];
    }
    return s;
}

void Vector::DangerousBuild(double *ptr, size_t s)
{
    t=ptr;
    nlin()=s;
    count=new int[1];
    *count=1;
}

void Vector::DangerousKill()
{
    t=0;
    nlin()=0;
    delete[] count;
}

Vector Vector::conv(const Vector& v) const {
    if (v.nlin()<nlin()) return v.conv(*this);

    Vector p(nlin()+v.nlin()-1);
    p.set(0);
    for (size_t i=0; i<v.nlin(); i++) {
#ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)((int)nlin(), v(i), t, 1, p.t+i, 1);
#else
        for (size_t j=0;j<nlin();j++)
        {
            p(i+j)+=v(i)*t[j];
        }
#endif
    }
    return p;
}

Vector Vector::conv_trunc(const Vector& v) const {
    Vector p(v.nlin());
    p.set(0);
    for (size_t i=0; i<v.nlin(); i++)
    {
        size_t m = std::min(nlin(),v.nlin()-i);
#ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)((int)m, v(i), t, 1, p.t+i, 1);
#else
        for (size_t j=0;j<m;j++)
        {
            p(i+j)+=v(i)*t[j];
        }
#endif
    }
    return p;
}


Matrix Vector::outer_product(const Vector& v) const {
    assert(nlin()==v.size());
    Matrix A(nlin(),nlin());
#ifdef HAVE_BLAS
    DGEMM(CblasNoTrans,CblasNoTrans,
        (int)nlin(),(int)nlin(),1,
        1.,t,(int)nlin(),
        v.t,(int)nlin(),
        0.,A.data(),(int)nlin());
#else
    for( unsigned int j=0; j<nlin(); j++ ) {
        for ( unsigned int i=0; i<nlin(); i++) {
            A(i,j)=v(i)*(*this)(j);
        }
    }
#endif
    return A;
}

void Vector::info() const {
    if (size() == 0) {
        std::cout << "Vector Empty" << std::endl;
        return;
    }

    std::cout << "Size : " << nlin() << std::endl;

    double minv = this->operator()(0);
    double maxv = this->operator()(0);
    size_t mini = 0;
    size_t maxi = 0;

    for(size_t i = 0; i < nlin(); ++i)
    {
        if (minv > this->operator()(i)) {
            minv = this->operator()(i);
            mini = i;
        } else if (maxv < this->operator()(i)) {
            maxv = this->operator()(i);
            maxi = i;
        }
    }

    std::cout << "Min Value : " << minv << " (" << mini << ")" << std::endl;
    std::cout << "Max Value : " << maxv << " (" << maxi << ")" << std::endl;
    std::cout << "First Values" << std::endl;
    for(size_t i = 0; i < std::min(nlin(),(size_t) 5); ++i) {
        std::cout << this->operator()(i) << std::endl;
    }
}

// =======
// = IOs =
// =======

void Vector::saveTxt(const char* filename) const
{
    maths::ofstream ofs(filename);
    ofs << maths::format("ascii") << *this;
}

void Vector::saveBin(const char* filename) const
{
    maths::ofstream ofs(filename);
    ofs << maths::format("old_binary") << *this;
}

void Vector::loadTxt(const char* filename)
{
    maths::ifstream ifs(filename);
    ifs >> maths::format("ascii") >> *this;
}

void Vector::loadBin(const char* filename)
{
    maths::ifstream ifs(filename);
    ifs >> maths::format("old_binary") >> *this;
}

void Vector::load( const char *filename ) {
    try {
        maths::ifstream ifs(filename);
        ifs >> *this;
    }
    catch (std::string s) {
        std::cout << s << std::endl;
    }
}

void Vector::save( const char *filename ) const {
    try {
        maths::ofstream ofs(filename);
        ofs << *this;
    }
    catch (std::string s) {
        std::cout << s << std::endl;
    }
}
