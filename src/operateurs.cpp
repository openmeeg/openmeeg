#include "operateurs.h"
#include "integrateur.h"
#include "analytiques.h"

//#define USE_SYM

#ifndef OPTIMIZED_OPERATOR_D
void operateurD(Geometry &geo,int I,int J,symmatrice &mat,int offsetI,int offsetJ)
// This function (NON OPTIMIZED VERSION) has the following arguments:
//    One geometry
//    the indices of the treated layers I and J
//    the storage matrix for the result
//    the upper left corner of the submatrix to be written is the matrix
{
    std::cout<<"OPERATEUR D..."<<std::endl;

    const Mesh &m1=geo.getM(I);
    const Mesh &m2=geo.getM(J);

    for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
        for(int j=offsetJ;j<offsetJ+m2.nbPts();j++)
        {
            // P1 functions are tested thus looping on vertices
            mat(i,j)=_operateurD(i-offsetI,j-offsetJ,m1,m2);
        }
}

void operateurD(Mesh &m1, Mesh &m2,matrice &mat,int offsetI,int offsetJ)
// This function (NON OPTIMIZED VERSION) has the following arguments:
//    One geometry
//    the indices of the treated layers I and J
//    the storage matrix for the result
//    the upper left corner of the submatrix to be written is the matrix
{
    std::cout<<"OPERATEUR D... (arg : mesh m1, mesh m2)"<<std::endl;

    for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
        for(int j=offsetJ;j<offsetJ+m2.nbPts();j++)
        {
            // P1 functions are tested thus looping on vertices
            mat(i,j)=_operateurD(i-offsetI,j-offsetJ,m1,m2);
        }
}

#else
void operateurD(Geometry &geo,int I,int J,symmatrice &mat,int offsetI,int offsetJ)
{
    // This function (OPTIMIZED VERSION) has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix

    std::cout<<"OPERATEUR D (Optimized)..."<<std::endl;

    const Mesh &m1=geo.getM(I);
    const Mesh &m2=geo.getM(J);

    for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
        for(int j=offsetJ;j<offsetJ+m2.nbTrgs();j++)
        {
            //In this version of the funtcion, in order to skip multiple computations of the same quantities
            //    loops are run over the triangles but the matrix cannot be filled in this function anymore
            //    That's why the filling is done is function _operateurD
            _operateurD(i-offsetI,j-offsetJ,m1,m2,mat,offsetI,offsetJ);
        }
}

void operateurD(Mesh &m1, Mesh &m2,genericMatrix &mat,int offsetI,int offsetJ)
{
    // This function (OPTIMIZED VERSION) has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix

    std::cout<<"OPERATEUR D (Optimized) ... (arg : mesh m1, mesh m2)"<<std::endl;

    for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
        for(int j=offsetJ;j<offsetJ+m2.nbTrgs();j++)
        {
            //In this version of the funtcion, in order to skip multiple computations of the same quantities
            //    loops are run over the triangles but the matrix cannot be filled in this function anymore
            //    That's why the filling is done is function _operateurD
            _operateurD(i-offsetI,j-offsetJ,m1,m2,mat,offsetI,offsetJ);
        }
}
#endif


void operateurS(Geometry &geo,int I,int J,symmatrice &mat,int offsetI,int offsetJ )
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

#ifndef USE_SYM
    if(I==J) {
        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbTrgs();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
        }
    } else {
        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
            // #pragma omp parallel for
            for(int j=offsetJ;j<offsetJ+m2.nbTrgs();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
        }
    }
#else
    if(I==J)
        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbTrgs();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
    else
    {
        if(m1.nbTrgs()<m2.nbTrgs())
            for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
                // #pragma omp parallel for
                for(int j=offsetJ+i-offsetI;j<offsetJ+m2.nbTrgs();j++)
                {
                    mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
        else
            for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
                // #pragma omp parallel for
                for(int j=offsetJ;j<=offsetJ+i-offsetI;j++)
                {
                    mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
    }
#endif
}

void operateurS(const Mesh &m1,const Mesh &m2,genericMatrix &mat,int offsetI,int offsetJ )
{
    // This function has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix

    std::cout<<"OPERATEUR S... (arg : mesh m1, mesh m2)"<<std::endl;

    // The operator S is given by Sij=\Int G*PSI(I,i)*Psi(J,j) with PSI(A,a) is a P0 test function on layer A and triangle a

#ifndef USE_SYM
    if(&m1==&m2)
        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbTrgs();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
    else
    {
        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
            for(int j=offsetJ;j<offsetJ+m2.nbTrgs();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
    }
#else

    if(I==J)
        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbTrgs();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
    else
    {
        if(m1.nbTrgs()<m2.nbTrgs())
            for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
                for(int j=offsetJ+i-offsetI;j<offsetJ+m2.nbTrgs();j++)
                {
                    mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
        else
            for(int i=offsetI;i<offsetI+m1.nbTrgs();i++)
                for(int j=offsetJ;j<=offsetJ+i-offsetI;j++)
                {
                    mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
    }
#endif
}

void operateurN(Geometry &geo,int I,int J,symmatrice &mat,int offsetI,int offsetJ , int IopS, int JopS )
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

#ifndef USE_SYM

    if(I==J)
        for(int i=offsetI;i<offsetI+m1.nbPts();i++)
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbPts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
    else
    {
        for(int i=offsetI;i<offsetI+m1.nbPts();i++)
            // #pragma omp parallel for
            for(int j=offsetJ;j<offsetJ+m2.nbPts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
    }

#else

    if(I==J)
        for(int i=offsetI;i<offsetI+m1.nbPts();i++)
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbPts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
    else
    {
        if(m1.nbPts()<m2.nbPts())
            for(int i=offsetI;i<offsetI+m1.nbPts();i++)
                // #pragma omp parallel for
                for(int j=offsetJ+i-offsetI;j<offsetJ+m2.nbPts();j++)
                {
                    mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
        else
            for(int i=offsetI;i<offsetI+m1.nbPts();i++)
                // #pragma omp parallel for
                for(int j=offsetJ;j<=offsetJ+i-offsetI;j++)
                {
                    mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
    }

#endif
}

void operateurN(Mesh &m1, Mesh &m2,genericMatrix &mat,int offsetI,int offsetJ,int IopS, int JopS)
{
    // This function has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix
    //  the upper left corner of the corresponding S block

    std::cout<<"OPERATEUR N... (arg : mesh m1, mesh m2)"<<std::endl;

#ifndef USE_SYM

    if(&m1==&m2)
        for(int i=offsetI;i<offsetI+m1.nbPts();i++)
            // #pragma omp parallel for schedule(runtime) ordered
            for(int j=i;j<offsetJ+m2.nbPts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
    else
    {
        for(int i=offsetI;i<offsetI+m1.nbPts();i++){
            for(int j=offsetJ;j<offsetJ+m2.nbPts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
        }
    }

#else

    if(I==J)
        for(int i=offsetI;i<offsetI+m1.nbPts();i++)
            for(int j=i;j<offsetJ+m2.nbPts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
    else
    {
        if(m1.nbPts()<m2.nbPts())
            for(int i=offsetI;i<offsetI+m1.nbPts();i++)
                for(int j=offsetJ+i-offsetI;j<offsetJ+m2.nbPts();j++)
                {
                    mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
        else
            for(int i=offsetI;i<offsetI+m1.nbPts();i++)
                for(int j=offsetJ;j<=offsetJ+i-offsetI;j++)
                {
                    mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
    }

#endif
}


//general routine for applying _operateurFerguson (see this function for further comments) to an entire mesh, and storing coordinates of the output in a matrix.
void operateurFerguson( const Vect3 x, const Mesh &m1, matrice &mat, int offsetI, int offsetJ)
{
    for(int j=offsetJ;j<offsetJ+m1.nbPts();j++)
    {
        Vect3 &v=_operateurFerguson(x,j-offsetJ,m1);
        mat(offsetI+0,j) = v.x();
        mat(offsetI+1,j) = v.y();
        mat(offsetI+2,j) = v.z();
    }
}

void operateurDipolePotDer(const Vect3 &r0, const Vect3 &q,const Mesh &inner_layer, vecteur &rhs, int offsetIdx)
{
    static analytiqueDipPotDer anaDPD;
    static integrateur<Vect3> gauss;
    gauss.setOrdre(GaussOrder);

    for(int i=0;i<inner_layer.nbTrgs();i++)
    {
        anaDPD.init(inner_layer,i,q,r0);
        Vect3 v=gauss.integre(anaDPD,inner_layer.getTrg(i),inner_layer);
        rhs(inner_layer.getTrg(i-offsetIdx).s1()+offsetIdx)+=v[0];
        rhs(inner_layer.getTrg(i-offsetIdx).s2()+offsetIdx)+=v[1];
        rhs(inner_layer.getTrg(i-offsetIdx).s3()+offsetIdx)+=v[2];
    }
}

void operateurDipolePot(const Vect3 &r0, const Vect3 &q, const Mesh &inner_layer, vecteur &rhs, int offsetIdx)
{
    static analytiqueDipPot anaDP;
    static integrateur<double> gauss;
    gauss.setOrdre(GaussOrder);

    anaDP.init(q,r0);
    for(int i=offsetIdx;i<offsetIdx+inner_layer.nbTrgs();i++)
    {
        rhs(i)+=gauss.integre(anaDP,inner_layer.getTrg(i-offsetIdx),inner_layer);
    }
}

// Grad wrt r0
void operateurDipolePotDerGrad(const Vect3 &r0, const Vect3 &q,const Mesh &inner_layer, vecteur rhs[6], int offsetIdx)
{
    static analytiqueDipPotDerGrad anaDPD;
    static integrateur< vect3array<6> > gauss;
    gauss.setOrdre(GaussOrder);

    for(int i=0;i<inner_layer.nbTrgs();i++)
    {
        anaDPD.init(inner_layer,i,q,r0);
        vect3array<6> v=gauss.integre(anaDPD,inner_layer.getTrg(i),inner_layer);
		for (int d=0;d<6;d++) { // derivatives wrt r0,q
			rhs[d](inner_layer.getTrg(i-offsetIdx).s1()+offsetIdx)+=v(d)[0];
			rhs[d](inner_layer.getTrg(i-offsetIdx).s2()+offsetIdx)+=v(d)[1];
			rhs[d](inner_layer.getTrg(i-offsetIdx).s3()+offsetIdx)+=v(d)[2];
		}
    }
}

// Grad wrt r0,q
void operateurDipolePotGrad(const Vect3 &r0, const Vect3 &q, const Mesh &inner_layer, vecteur rhs[6], int offsetIdx)
{
    static analytiqueDipPotGrad anaDP;
    static integrateur< vect3array<2> > gauss;
    gauss.setOrdre(GaussOrder);

    anaDP.init(q,r0);
    for(int i=offsetIdx;i<offsetIdx+inner_layer.nbTrgs();i++)
    {
        vect3array<2> v = gauss.integre(anaDP,inner_layer.getTrg(i-offsetIdx),inner_layer);
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
