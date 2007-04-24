/*! \file
\brief File defining the template TPoint
*/
#ifndef _TPoint_H
#define _TPoint_H

#include "utils.h"

/*! \addtogroup TPointGroup TPoint and associated functions
    \brief Here is a list of functions close to TPoint.
*/
/*! Template for points
    \ingroup TPointGroup
    \param T type of each coordinate of the point (\c int / \c double / \c reel / ... )
    \param dim dimension of the point
*/

/*! \return the product of two \c reel.
    \ingroup reelGroup
*/
inline reel mult( const reel &a, const reel &b) {return a*b;}

template <typename T,int dim> class TPoint ;

/*! multiply two TPoint.
    \ingroup TPointGroup
*/
template <typename T,int dim> inline TPoint<T,dim> mult(const TPoint<T,dim> &x, const TPoint<T,dim>& p) { // coord a coord.
    return x.mult(p);
}

/*! divide two TPoint
    \ingroup TPointGroup
*/
template <typename T,int dim> inline TPoint<T,dim> div(const TPoint<T,dim> &x, const TPoint<T,dim>& p) { // coord a coord.
    return x.div(p);
}

/*! \return the division of two \c reel.
    \ingroup reelGroup
*/
inline reel div( const reel &a, const reel &b) { return a/b; }

template <typename T,int dim> class TPoint
{
protected:
    T x[dim]; //!< Array containing the data of the point
public:
    //! Constructor.
    inline TPoint(){}
    //! Copy constructor.
    inline TPoint(const T X[dim]){
        int i;
        for (i=0;i<dim;i++)
            x[i]=X[i];
    }
    //! Copy contructor
    template <typename Y> inline TPoint(const TPoint<Y,dim>& c)
    {
        for (int i=0;i<dim;i++)
            x[i]=T(c[i]);
    }

    /*! \brief Fast initialisation of a TPoint.

        Usually used to initialise a TPoint to 0.
        \code
        //declare P
        TPoint<double, 3> P;
        //initialise P
        P.set(0);
        \endcode
    */
    inline void set(T X) {
        int i;
        for (i=0;i<dim;i++)
            x[i]=X;
    }

    //! Fast initialisation of a 2D point.
    inline TPoint(T a,T b) { assert(dim==2) ; x[0]=a;x[1]=b; }
    //! Fast initialisation of a 3D point.
    inline TPoint(T a,T b,T c) { assert(dim==3); x[0]=a;x[1]=b;x[2]=c; }
    //! Copy of a point.
    inline TPoint(T X) { set(X); }

    /*! Set the value of a coordinate.
        \param d index of the coordinate
        \param x new value
    */
    inline TPoint<T,dim> subst(int d,T xx) const {
        TPoint<T,dim> p(*this);
        p[d]=xx;
        return p;
    }

    /*! Set the value of coordinates [<em>0 ... d</em>]
        \param d index of the final coordinate.
        \param X TPoint containing the new values.
    */
    inline TPoint<T,dim> subst(int d,TPoint<T,dim> X) const {
        TPoint<T,dim> p;
        for (int i=0;i<dim;i++)
            p[i]=(i<=d)?X[i]:x[i];
        return p;
    }

    //! Operator[] (read)
    inline T operator[](int i) const {
        assert(i<dim);
        return x[i];
    }

    //! Operator[] (write)
    inline T& operator[](int i) {
        assert(i<dim);
        return x[i];
    }

    //! Equality test
    inline bool operator==(const TPoint<T,dim>& p) const {
        int i;
        for (i=0;i<dim;i++)
            if (x[i]!=p[i])
                return false;
        return true;
    }

    //! Difference test
    inline bool operator!=(const TPoint<T,dim>& p) const {
        return !(*this == p);
    }


    //! Copy the values of a TPoint to another.
    inline void operator=(const TPoint<T,dim>& p){
        for (int i=0;i<dim;i++)
            x[i]=p[i];
    }

    //! Add two TPoints coordinate by coordinate.
    inline TPoint<T,dim> operator+(const TPoint<T,dim>& p) const {
        TPoint<T,dim> s;
        for (int i=0;i<dim;i++)
            s[i]=x[i]+p[i];
        return s;
    }

    //! Add a TPoint to a TPoint coordinate by coordinate.
    inline void operator+=(const TPoint<T,dim>& p){
        for (int i=0;i<dim;i++)
            x[i] += p[i];
    }

    //! Substract two TPoints coordinate by coordinate.
    inline TPoint<T,dim> operator-(const TPoint<T,dim>& p) const {
        TPoint<T,dim> s;
        for (int i=0;i<dim;i++)
            s[i]=x[i]-p[i];
        return s;
    }

    //! Multiply the TPoint by -1.
    inline TPoint<T,dim> operator-() const {
        TPoint<T,dim> s;
        for (int i=0;i<dim;i++)
            s[i]=-x[i];
        return s;
    }

    /*! Compute a^b (dimension 2).
        \return a[0]*b[1]-a[1]*b[0] . The return value is a \c reel.
    */
    inline T operator^(const TPoint<T,2>& p) const {
        assert(dim==2);
        return x[0]*p[1]-x[1]*p[0];
    }

    /*! \brief Compute a^b in dimension 1 (ERROR)
        Useless for user, useful to remove compiling warning/error.
        \warning This method should never be used !
    */
    inline T operator^(const TPoint<T,1>& p) const {
        assert(false); // This method should never be used
        return p[0];
    }

    /*! Compute the vectorial product of two TPoint<T,3>.
        \attention dimension of the TPoint has to be 3.
        \return s=a^b , dimension of \c s is 3.
    */
    inline TPoint<T,3> operator^(const TPoint<T,3>& p) const {
        TPoint<T,3> s;
        s[0]=x[1]*p[2]-x[2]*p[1];
        s[1]=x[2]*p[0]-x[0]*p[2];
        s[2]=x[0]*p[1]-x[1]*p[0];
        return s;
    }

    /*! Compute the scalar product of two TPoint.
        \return the return value is a \c reel .
    */
    inline T operator*(const TPoint<T,dim>& p) const {
        T s=0.0;
        for (int i=0;i<dim;i++)
            s=s+::mult(x[i],p[i]);
        return s;
    }

    /*! Compute the product of a TPoint by a scalar.
        \param v scalar.
        \return point * v.
    */
    inline TPoint<T,dim> operator*(T v) const {
        TPoint<T,dim> s;
        for (int i=0;i<dim;i++)
            s[i]=(x[i]*v);
        return s;
    }

    /*! Compute the division of a TPoint by a scalar.
        \param v scalar.
        \return point / v.
    */
    inline TPoint<T,dim> operator/(T v) const {
        TPoint<T,dim> s;
        for (int i=0;i<dim;i++)
            s[i]=(x[i]/v);
        return s;
    }

    /*! Divide a TPoint by a scalar for each coordinate.
        \param v scalar
        \return void
        \attention v should be non-zero...
    */
    inline void operator/=(T v) {
        for (int i=0;i<dim;i++)
            x[i] /= v;
    }

    /*! Compute the product of a TPoint by a TPoint for each coordinate.
        \param p Tpoint.
        \return s = x.mult(p) <=> s[0]=x[0]*p[0] . s[1]=x[1]*p[1]. etc.
    */
    inline TPoint<T,dim> mult(const TPoint<T,dim>& p) const {
        TPoint<T,dim> s;
        for (int i=0;i<dim;i++)
            s[i]=::mult(x[i],p[i]);
        return s;
    }

    /*! Compute the division of a TPoint by a TPoint for each coordinate.
        \param p Tpoint.
        \return s = x.div(p) <=> s[0]=x[0]/p[0] . s[1]=x[1]/p[1]. etc.
    */
    inline TPoint<T,dim> div(const TPoint<T,dim>& p) const {
        TPoint<T,dim> s;
        for (int i=0;i<dim;i++)
            s[i]=::div(x[i],p[i]);
        return s;
    }

    /*! Compute the product of each coordinate of a \c TPoint<T,dim> by a scalar whose type is not \c T .
    */
    template <typename S> inline TPoint<T,dim> times(S v) const { // product by a scalar whose type is not \c T
                                                                  // e.g: reel x;
                                                                  //      TPoint<reel,dim> p,q;
                                                                  //      q=p*x; // OK
                                                                  //      TPoint<dTPoint,dim> P,Q;
                                                                  //      Q=P*x; // bad: convert x to dTPoint!
                                                                  //      Q=P.times(x); // OK
        TPoint<T,dim> s;
        for (int i=0;i<dim;i++)
            s[i]=T(x[i]*v);
        return s;
    }

    /*! \return the squared norm of a TPoint.
        \attention The type of return value is T.\n
        Beware to overflow with TPoint<char,dim> and TPoint<short,dim>.\n
        Prefer dnorme2().
    */
    inline T norme2() const {
        T n=0;
        for (int i=0;i<dim;i++)
            n+=x[i]*x[i];
        return n;
    }

    /*! \return the squared norm of a TPoint as a \c int .
        \attention Less sensible to overflow than norme2().
    */
    inline int inorme2() const {
        int n=0;
        for (int i=0;i<dim;i++)
            n+=int(x[i])*x[i];
        return n;
    }

    /*! \return the squared norm of a TPoint as a \c reel .
        \attention Less sensible to overflow than norme2().
    */
    inline reel dnorme2() const {
        reel n=0;
        for (int i=0;i<dim;i++)
            n+=reel(x[i])*x[i];
        return n;
    }

    /*! \return the norm of a \c TPoint. Return type is \c reel.
    */
    inline reel norme() const {
        return (sqrt(dnorme2()));
    }

    /*! Normalize a vector.
        \attention If the TPoint == 0, it doesn't change.
    */
    inline void normalize() {
        reel n=norme();
        if (n==0)
            return;
        for (int i=0;i<dim;i++)
            x[i]=T(x[i]/n);
    }

    /*! Normalized version of a vector.
    */
    inline TPoint<T,dim> normalized() const {
        reel n=norme();
        if (n==0)
            return *this;
        TPoint<T,dim> p;
        for (int i=0;i<dim;i++)
            p.x[i]=T(x[i]/n);
        return p;
    }

    /*! \return the sum of coordinates
    */
    inline T sum() const {
        T n=0;
        for (int i=0;i<dim;i++)
            n+=x[i];
        return n;
    }

    /*! \return the product of coordinates
        \attention Since the return type is T, it may be some overflows.\n
        If possible prefer coords::selfprod().
    */
    inline T prod() const {
        T p;
        int i;
        for (i=0,p=1;i<dim;i++)
            p*=x[i];
        return p;
    }

    //! Binary write of a TPoint to a stream.
    inline void write(std::ostream& f) const {
        f.write((char*)this->x,dim*sizeof(T));
    }

    //! Binary read of a TPoint from a stream.
    inline void read(std::istream& f) {
        f.read((char*)this->x,dim*sizeof(T));
    }

    //! Ascii write of a TPoint to a stream.
    friend inline std::ostream& operator<<(std::ostream& f,const TPoint<T,dim>& p) {
        for (int i=0;i<dim;i++)
            f<<p.x[i]<<((i<dim-1)?" ":"");
        return f;
    }

    //! Ascii read of a TPoint from a stream.
    friend inline std::istream& operator>>(std::istream& f,TPoint<T,dim>& p) {
        for (int i=0;i<dim;i++)
            f>>p.x[i];
        return f;
    }

    //! Left-multiply a TPoint by a scalar.
    friend inline TPoint<T,dim> operator*(T v,const TPoint<T,dim>& p) {
        return p*v;
    }

};

/*! \return the absolute value of each coordinate of a TPoint.
    \ingroup TPointGroup
*/
template <typename T,int dim> inline TPoint<T,dim> abs(const TPoint<T,dim>& m) {
    TPoint<T,dim> s;
    for (int i=0;i<dim;i++)
        s[i]=std::max(m[i],-m[i]);
    return s;
}

/*! \return the maximum of each coordinate of two TPoint.
    \ingroup TPointGroup
*/
template <typename T,int dim> inline TPoint<T,dim> maxi(const TPoint<T,dim>& m,const TPoint<T,dim> &n){
    TPoint<T,dim> s;
    for (int i=0;i<dim;i++)
        s[i]=std::max(m[i],n[i]);
    return s;
}

/*! \return the maximum of the coordinates of a TPoint.
    \ingroup TPointGroup
*/
template <typename T,int dim> inline T maxi(const TPoint<T,dim>& p){
    T m=p[0];
    for (int i=1;i<dim;i++)
        m=std::max(m,p[i]);
    return m;
}

/*! \return the minimum of each coordinate of two TPoint.
    \ingroup TPointGroup
*/
template <typename T,int dim> inline TPoint<T,dim> mini(const TPoint<T,dim>& m,const TPoint<T,dim> &n){
    TPoint<T,dim> s;
    for (int i=0;i<dim;i++)
        s[i]=std::min(m[i],n[i]);
    return s;
}

/*! \return the minimum of the coordinates of a TPoint.
    \ingroup TPointGroup
*/
template <typename T,int dim> inline T mini(const TPoint<T,dim>& p){
    T m=p[0];
    for (int i=1;i<dim;i++)
        m=std::min(m,p[i]);
    return m;
}

/*! rPoint(dim) is a faster way to write TPoint<reel,dim> , a TPoint with coordinates of type \c reel.
    \param dim dimension of the point
    \ingroup TPointGroup
*/
#define rPoint(dim) TPoint<reel,dim>

/*! iPoint(dim) is a faster way to write TPoint<int,dim> , a TPoint with coordinates of type \c int.
    \param dim dimension of the point
    \ingroup TPointGroup
*/
#define iPoint(dim) TPoint<int,dim>

/*! sPoint(dim) is a faster way to write TPoint<short,dim> , a TPoint with coordinates of type \c short.
    \param dim dimension of the point
    \ingroup TPointGroup
*/
#define sPoint(dim) TPoint<short,dim>

#endif

