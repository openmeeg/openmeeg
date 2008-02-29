/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#include "operators.h"
#include "integrator.h"
#include "analytics.h"

#ifndef OPTIMIZED_OPERATOR_D
void operatorD(const Geometry &geo,const int I,const int J,const int GaussOrder,symmatrice &mat,const int offsetI,const int offsetJ)
// This function (NON OPTIMIZED VERSION) has the following arguments:
//    One geometry
//    the indices of the treated layers I and J
//    the storage matrix for the result
//    the upper left corner of the submatrix to be written is the matrix
{
    std::cout<<"OPERATEUR D..."<<std::endl;

    const Mesh &m1=geo.getM(I);
    const Mesh &m2=geo.getM(J);

    for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
        progressbar(i-offsetI,m1.nbTrgs());
        #ifdef USE_OMP
        #pragma omp parallel for
        #endif
        for(int j=offsetJ;j<offsetJ+m2.nbPts();j++)
        {
            // P1 functions are tested thus looping on vertices
            mat(i,j)=_operatorD(i-offsetI,j-offsetJ,GaussOrder,m1,m2);
        }
    }
}

void operatorD(const Mesh &m1,const Mesh &m2,genericMatrix &mat,const int offsetI,const int offsetJ,const int GaussOrder)
// This function (NON OPTIMIZED VERSION) has the following arguments:
//    One geometry
//    the indices of the treated layers I and J
//    the storage matrix for the result
//    the upper left corner of the submatrix to be written is the matrix
{
    std::cout<<"OPERATEUR D... (arg : mesh m1, mesh m2)"<<std::endl;

    for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
        progressbar(i-offsetI,m1.nbTrgs());
        #ifdef USE_OMP
        #pragma omp parallel for
        #endif
        for(int j=offsetJ;j<offsetJ+m2.nbPts();j++)
        {
            // P1 functions are tested thus looping on vertices
            mat(i,j)=_operatorD(i-offsetI,j-offsetJ,GaussOrder,m1,m2);
        }
    }
}

#else // OPTIMIZED_OPERATOR_D

void operatorD(const Geometry &geo,const int I,const int J,const int GaussOrder,symmatrice &mat,const int offsetI,const int offsetJ)
{
    // This function (OPTIMIZED VERSION) has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix

    std::cout<<"OPERATEUR D (Optimized)..."<<std::endl;

    const Mesh &m1=geo.getM(I);
    const Mesh &m2=geo.getM(J);

    for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
        progressbar(i-offsetI,m1.nbTrgs());
        for(int j=offsetJ;j<offsetJ+m2.nbTrgs();j++)
        {
            //In this version of the funtcion, in order to skip multiple computations of the same quantities
            //    loops are run over the triangles but the matrix cannot be filled in this function anymore
            //    That's why the filling is done is function _operatorD
            _operatorD(i-offsetI,j-offsetJ,GaussOrder,m1,m2,mat,offsetI,offsetJ);
        }
    }
}

void operatorD(const Mesh &m1,const Mesh &m2,genericMatrix &mat,const int offsetI,const int offsetJ,const int GaussOrder)
{
    // This function (OPTIMIZED VERSION) has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix

    std::cout<<"OPERATEUR D (Optimized) ... (arg : mesh m1, mesh m2)"<<std::endl;

    for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
        progressbar(i-offsetI,m1.nbTrgs());
        for(int j=offsetJ;j<offsetJ+m2.nbTrgs();j++)
        {
            //In this version of the funtcion, in order to skip multiple computations of the same quantities
            //    loops are run over the triangles but the matrix cannot be filled in this function anymore
            //    That's why the filling is done is function _operatorD
            _operatorD(i-offsetI,j-offsetJ,GaussOrder,m1,m2,mat,offsetI,offsetJ);
        }
    }
}

#endif // OPTIMIZED_OPERATOR_D

void operatorS(const Geometry &geo,const int I,const int J,const int GaussOrder,symmatrice &mat,const int offsetI,const int offsetJ )
{
    // This function has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix

    std::cout<<"OPERATEUR S..."<<std::endl;

    // The operator S is given by Sij=\Int G*PSI(I,i)*Psi(J,j) with PSI(A,a) is a P0 test function on layer A and triangle a
    const Mesh &m1=geo.getM(I);
    const Mesh &m2=geo.getM(J);

    if(I==J) {
        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
            progressbar(i-offsetI,m1.nbTrgs());
            #ifdef USE_OMP
            #pragma omp parallel for
            #endif
            for(int j=i;j<offsetJ+m2.nbTrgs();j++)
            {
                mat(i,j)=_operatorS(i-offsetI,j-offsetJ,GaussOrder,m1,m2);
            }
        }
    } else {
        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
            progressbar(i-offsetI,m1.nbTrgs());
            #ifdef USE_OMP
            #pragma omp parallel for
            #endif
            for(int j=offsetJ;j<offsetJ+m2.nbTrgs();j++)
            {
                mat(i,j)=_operatorS(i-offsetI,j-offsetJ,GaussOrder,m1,m2);
            }
        }
    }
}

void operatorS(const Mesh &m1,const Mesh &m2,genericMatrix &mat,const int offsetI,const int offsetJ,const int GaussOrder)
{
    // This function has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix

    std::cout<<"OPERATEUR S... (arg : mesh m1, mesh m2)"<<std::endl;

    // The operator S is given by Sij=\Int G*PSI(I,i)*Psi(J,j) with PSI(A,a) is a P0 test function on layer A and triangle a
    if(&m1==&m2)
        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
            progressbar(i-offsetI,m1.nbTrgs());
            #ifdef USE_OMP
            #pragma omp parallel for
            #endif
            for(int j=i;j<offsetJ+m2.nbTrgs();j++)
            {
                mat(i,j)=_operatorS(i-offsetI,j-offsetJ,GaussOrder,m1,m2);
            }
        }
    else
    {
        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
            progressbar(i-offsetI,m1.nbTrgs());
            #ifdef USE_OMP
            #pragma omp parallel for
            #endif
            for(int j=offsetJ;j<offsetJ+m2.nbTrgs();j++)
            {
                mat(i,j)=_operatorS(i-offsetI,j-offsetJ,GaussOrder,m1,m2);
            }
        }
    }
}

void operatorN(const Geometry &geo,const int I,const int J,const int GaussOrder,symmatrice &mat,const int offsetI,const int offsetJ,const int IopS,const int JopS)
{
    // This function has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix
    //  the upper left corner of the corresponding S block
    std::cout<<"OPERATEUR N..."<<std::endl;

    const Mesh &m1=geo.getM(I);
    const Mesh &m2=geo.getM(J);

    if(I==J)
        for(int i=offsetI;i<offsetI+m1.nbPts();i++) {
            progressbar(i-offsetI,m1.nbPts());
            #ifdef USE_OMP
            #pragma omp parallel for
            #endif
            for(int j=i;j<offsetJ+m2.nbPts();j++)
            {
                mat(i,j)=_operatorN(i-offsetI,j-offsetJ,GaussOrder,m1,m2,IopS,JopS,mat);
            }
        }
    else
    {
        for(int i=offsetI;i<offsetI+m1.nbPts();i++) {
            progressbar(i-offsetI,m1.nbPts());
            #ifdef USE_OMP
            #pragma omp parallel for
            #endif
            for(int j=offsetJ;j<offsetJ+m2.nbPts();j++)
            {
                mat(i,j)=_operatorN(i-offsetI,j-offsetJ,GaussOrder,m1,m2,IopS,JopS,mat);
            }
        }
    }
}

void operatorN(const Mesh &m1,const Mesh &m2,genericMatrix &mat,const int offsetI,const int offsetJ,const int GaussOrder,const int IopS,const int JopS)
{
    // This function has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix
    //  the upper left corner of the corresponding S block

    std::cout<<"OPERATEUR N... (arg : mesh m1, mesh m2)"<<std::endl;

    if(&m1==&m2) {
        for(int i=offsetI;i<offsetI+m1.nbPts();i++) {
            progressbar(i-offsetI,m1.nbPts());
            #ifdef USE_OMP
            #pragma omp parallel for
            #endif
            for(int j=i;j<offsetJ+m2.nbPts();j++)
            {
                mat(i,j)=_operatorN(i-offsetI,j-offsetJ,GaussOrder,m1,m2,IopS,JopS,mat);
            }
        }
    } else {
        for(int i=offsetI;i<offsetI+m1.nbPts();i++){
            progressbar(i-offsetI,m1.nbPts());
            #ifdef USE_OMP
            #pragma omp parallel for
            #endif
            for(int j=offsetJ;j<offsetJ+m2.nbPts();j++)
            {
                mat(i,j)=_operatorN(i-offsetI,j-offsetJ,GaussOrder,m1,m2,IopS,JopS,mat);
            }
        }
    }
}

// general routine for applying _operatorFerguson (see this function for further comments)
// to an entire mesh, and storing coordinates of the output in a matrix.
void operatorFerguson(const Vect3& x, const Mesh &m1, matrice &mat, int offsetI, int offsetJ)
{
    #ifdef USE_OMP
    #pragma omp parallel for
    #endif
    for(int j=offsetJ;j<offsetJ+m1.nbPts();j++)
    {
        Vect3 v = _operatorFerguson(x,j-offsetJ,m1);
        mat(offsetI+0,j) = v.x();
        mat(offsetI+1,j) = v.y();
        mat(offsetI+2,j) = v.z();
    }
}

void operatorDipolePotDer(const Vect3 &r0,const Vect3 &q,const Mesh &inner_layer,vecteur &rhs,const int offsetIdx,const int GaussOrder)
{
    static analyticDipPotDer anaDPD;

#ifdef ADAPT_RHS
    adaptive_integrator<Vect3> gauss(0.001);
#else
    static integrator<Vect3> gauss;
#endif //ADAPT_RHS
    gauss.setOrder(GaussOrder);
    #ifdef USE_OMP
    #pragma omp parallel for private(anaDPD)
    #endif
    for(int i=0;i<inner_layer.nbTrgs();i++)
    {
        anaDPD.init(inner_layer,i,q,r0);
        Vect3 v=gauss.integrate(anaDPD,inner_layer.getTrg(i),inner_layer);
        #ifdef USE_OMP
        #pragma omp critical
        #endif
        {
        rhs(inner_layer.getTrg(i-offsetIdx).s1()+offsetIdx)+=v(0);
        rhs(inner_layer.getTrg(i-offsetIdx).s2()+offsetIdx)+=v(1);
        rhs(inner_layer.getTrg(i-offsetIdx).s3()+offsetIdx)+=v(2);
        }
    }
}

void operatorDipolePot(const Vect3 &r0, const Vect3 &q, const Mesh &inner_layer, vecteur &rhs,const int offsetIdx,const int GaussOrder)
{
    static analyticDipPot anaDP;

    anaDP.init(q,r0);
#ifdef ADAPT_RHS
    adaptive_integrator<double> gauss(0.001);
#else
    static integrator<double> gauss;
#endif
    gauss.setOrder(GaussOrder);
    #ifdef USE_OMP
    #pragma omp parallel for
    #endif
    for(int i=offsetIdx;i<offsetIdx+inner_layer.nbTrgs();i++)
    {
        double d = gauss.integrate(anaDP,inner_layer.getTrg(i-offsetIdx),inner_layer);
        #ifdef USE_OMP
        #pragma omp critical
        #endif
        rhs(i) += d;
    }
}

// Grad wrt r0
void operatorDipolePotDerGrad(const Vect3 &r0, const Vect3 &q,const Mesh &inner_layer, vecteur rhs[6],const int offsetIdx,const int GaussOrder)
{
    static analyticDipPotDerGrad anaDPD;
    static integrator< vect3array<6> > gauss;
    gauss.setOrder(GaussOrder);

    for(int i=0;i<inner_layer.nbTrgs();i++)
    {
        anaDPD.init(inner_layer,i,q,r0);
        vect3array<6> v=gauss.integrate(anaDPD,inner_layer.getTrg(i),inner_layer);
        for (int d=0;d<6;d++) { // derivatives wrt r0,q
            rhs[d](inner_layer.getTrg(i-offsetIdx).s1()+offsetIdx)+=v(d)(0);
            rhs[d](inner_layer.getTrg(i-offsetIdx).s2()+offsetIdx)+=v(d)(1);
            rhs[d](inner_layer.getTrg(i-offsetIdx).s3()+offsetIdx)+=v(d)(2);
        }
    }
}

// Grad wrt r0,q
void operatorDipolePotGrad(const Vect3 &r0,const Vect3 &q,const Mesh &inner_layer,vecteur rhs[6],const int offsetIdx,const int GaussOrder)
{
    static analyticDipPotGrad anaDP;
    static integrator< vect3array<2> > gauss;
    gauss.setOrder(GaussOrder);

    anaDP.init(q,r0);
    for(int i=offsetIdx;i<offsetIdx+inner_layer.nbTrgs();i++)
    {
        vect3array<2> v = gauss.integrate(anaDP,inner_layer.getTrg(i-offsetIdx),inner_layer);
        // grad_r0
        rhs[0](i) += v(0).x();
        rhs[1](i) += v(0).y();
        rhs[2](i) += v(0).z();
        // grad_q
        rhs[3](i) += v(1).z();
        rhs[4](i) += v(1).y();
        rhs[5](i) += v(1).z();
    }
}

void operatorP1P0(const Geometry &geo,const int I,symmatrice &mat,const int offsetI,const int offsetJ)
{
// This time mat(i,j)+= ... the matrix is incremented by the P1P0 operator
	std::cout<<"OPERATEUR P1P0..."<<std::endl;
	const Mesh &m=geo.getM(I);
	for(int i=offsetI;i<offsetI+m.nbTrgs();i++)
		for(int j=offsetJ;j<offsetJ+m.nbPts();j++)
		{
	   	mat(i,j)=-mat(i,j)+_operateurP1P0(i-offsetI,j-offsetJ,m)/2;
		}
}

