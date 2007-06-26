#include "operateurs.h"
#include "integrateur.h"
#include "analytiques.h"

//#define USE_SYM

#ifndef OPTIMIZED_OPERATOR_D
void operateurD(geometry &geo,int I,int J,symmatrice &mat,int offsetI,int offsetJ)
// This function (NON OPTIMIZED VERSION) has the following arguments:
//    One geometry
//    the indices of the treated layers I and J
//    the storage matrix for the result
//    the upper left corner of the submatrix to be written is the matrix
{
    std::cout<<"OPERATEUR D..."<<std::endl;

    const mesh &m1=geo.getM(I);
    const mesh &m2=geo.getM(J);

    for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
        for(int j=offsetJ;j<offsetJ+m2.nbr_pts();j++)
        {
            // P1 functions are tested thus looping on vertices
            mat(i,j)=_operateurD(i-offsetI,j-offsetJ,m1,m2);
        }
}

void operateurD(mesh &m1, mesh &m2,matrice &mat,int offsetI,int offsetJ)
// This function (NON OPTIMIZED VERSION) has the following arguments:
//    One geometry
//    the indices of the treated layers I and J
//    the storage matrix for the result
//    the upper left corner of the submatrix to be written is the matrix
{
    std::cout<<"OPERATEUR D... (arg : mesh m1, mesh m2)"<<std::endl;

    for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
        for(int j=offsetJ;j<offsetJ+m2.nbr_pts();j++)
        {
            // P1 functions are tested thus looping on vertices
            mat(i,j)=_operateurD(i-offsetI,j-offsetJ,m1,m2);
        }
}

#else
void operateurD(geometry &geo,int I,int J,symmatrice &mat,int offsetI,int offsetJ)
{
    // This function (OPTIMIZED VERSION) has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix

    std::cout<<"OPERATEUR D (Optimized)..."<<std::endl;

    const mesh &m1=geo.getM(I);
    const mesh &m2=geo.getM(J);

    for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
        for(int j=offsetJ;j<offsetJ+m2.nbr_trg();j++)
        {
            //In this version of the funtcion, in order to skip multiple computations of the same quantities
            //    loops are run over the triangles but the matrix cannot be filled in this function anymore
            //    That's why the filling is done is function _operateurD
            _operateurD(i-offsetI,j-offsetJ,m1,m2,mat,offsetI,offsetJ);
        }
}

void operateurD(mesh &m1, mesh &m2,genericMatrix &mat,int offsetI,int offsetJ)
{
    // This function (OPTIMIZED VERSION) has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix

    std::cout<<"OPERATEUR D (Optimized) ... (arg : mesh m1, mesh m2)"<<std::endl;

    for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
        for(int j=offsetJ;j<offsetJ+m2.nbr_trg();j++)
        {
            //In this version of the funtcion, in order to skip multiple computations of the same quantities
            //    loops are run over the triangles but the matrix cannot be filled in this function anymore
            //    That's why the filling is done is function _operateurD
            _operateurD(i-offsetI,j-offsetJ,m1,m2,mat,offsetI,offsetJ);
        }
}
#endif


void operateurS(geometry &geo,int I,int J,symmatrice &mat,int offsetI,int offsetJ )
{
    // This function has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix

    std::cout<<"OPERATEUR S..."<<std::endl;

    // The operator S is given by Sij=\Int G*PSI(I,i)*Psi(J,j) with PSI(A,a) is a P0 test function on layer A and triangle a
    const mesh &m1=geo.getM(I);
    const mesh &m2=geo.getM(J);

#ifndef USE_SYM
    if(I==J) {
        for(int i=offsetI;i<offsetI+m1.nbr_trg();i++) {
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbr_trg();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
        }
    } else {
        for(int i=offsetI;i<offsetI+m1.nbr_trg();i++) {
            // #pragma omp parallel for
            for(int j=offsetJ;j<offsetJ+m2.nbr_trg();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
        }
    }
#else
    if(I==J)
        for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbr_trg();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
    else
    {
        if(m1.nbr_trg()<m2.nbr_trg())
            for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
                // #pragma omp parallel for
                for(int j=offsetJ+i-offsetI;j<offsetJ+m2.nbr_trg();j++)
                {
                    mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
        else
            for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
                // #pragma omp parallel for
                for(int j=offsetJ;j<=offsetJ+i-offsetI;j++)
                {
                    mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
    }
#endif
}

void operateurS(const mesh &m1,const mesh &m2,genericMatrix &mat,int offsetI,int offsetJ )
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
        for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbr_trg();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
    else
    {
        for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
            for(int j=offsetJ;j<offsetJ+m2.nbr_trg();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
    }
#else

    if(I==J)
        for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbr_trg();j++)
            {
                mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
            }
    else
    {
        if(m1.nbr_trg()<m2.nbr_trg())
            for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
                for(int j=offsetJ+i-offsetI;j<offsetJ+m2.nbr_trg();j++)
                {
                    mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
        else
            for(int i=offsetI;i<offsetI+m1.nbr_trg();i++)
                for(int j=offsetJ;j<=offsetJ+i-offsetI;j++)
                {
                    mat(i,j)=_operateurS(i-offsetI,j-offsetJ,m1,m2);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
    }
#endif
}

void operateurN(geometry &geo,int I,int J,symmatrice &mat,int offsetI,int offsetJ , int IopS, int JopS )
{
    // This function has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage matrix for the result
    //    the upper left corner of the submatrix to be written is the matrix
    //  the upper left corner of the corresponding S block
    std::cout<<"OPERATEUR N..."<<std::endl;

    const mesh &m1=geo.getM(I);
    const mesh &m2=geo.getM(J);

#ifndef USE_SYM

    if(I==J)
        for(int i=offsetI;i<offsetI+m1.nbr_pts();i++)
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbr_pts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
    else
    {
        for(int i=offsetI;i<offsetI+m1.nbr_pts();i++)
            // #pragma omp parallel for
            for(int j=offsetJ;j<offsetJ+m2.nbr_pts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
    }

#else

    if(I==J)
        for(int i=offsetI;i<offsetI+m1.nbr_pts();i++)
            // #pragma omp parallel for
            for(int j=i;j<offsetJ+m2.nbr_pts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
    else
    {
        if(m1.nbr_pts()<m2.nbr_pts())
            for(int i=offsetI;i<offsetI+m1.nbr_pts();i++)
                // #pragma omp parallel for
                for(int j=offsetJ+i-offsetI;j<offsetJ+m2.nbr_pts();j++)
                {
                    mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
        else
            for(int i=offsetI;i<offsetI+m1.nbr_pts();i++)
                // #pragma omp parallel for
                for(int j=offsetJ;j<=offsetJ+i-offsetI;j++)
                {
                    mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
    }

#endif
}

void operateurN(mesh &m1, mesh &m2,genericMatrix &mat,int offsetI,int offsetJ,int IopS, int JopS)
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
        for(int i=offsetI;i<offsetI+m1.nbr_pts();i++)
            // #pragma omp parallel for schedule(runtime) ordered
            for(int j=i;j<offsetJ+m2.nbr_pts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
    else
    {
        for(int i=offsetI;i<offsetI+m1.nbr_pts();i++){
            for(int j=offsetJ;j<offsetJ+m2.nbr_pts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
        }
    }

#else

    if(I==J)
        for(int i=offsetI;i<offsetI+m1.nbr_pts();i++)
            for(int j=i;j<offsetJ+m2.nbr_pts();j++)
            {
                mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
            }
    else
    {
        if(m1.nbr_pts()<m2.nbr_pts())
            for(int i=offsetI;i<offsetI+m1.nbr_pts();i++)
                for(int j=offsetJ+i-offsetI;j<offsetJ+m2.nbr_pts();j++)
                {
                    mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
        else
            for(int i=offsetI;i<offsetI+m1.nbr_pts();i++)
                for(int j=offsetJ;j<=offsetJ+i-offsetI;j++)
                {
                    mat(i,j)=_operateurN(i-offsetI,j-offsetJ,m1,m2,IopS,JopS,mat);
                    mat(offsetI+j-offsetJ,offsetJ+i-offsetI)=mat(i,j);
                }
    }

#endif
}


//general routine for applying _operateurFerguson (see this function for further comments) to an entire mesh, and storing coordinates of the output in a matrix.
void operateurFerguson( const vect3 x, const mesh &m1, matrice &mat, int offsetI, int offsetJ)
{
    for(int j=offsetJ;j<offsetJ+m1.nbr_pts();j++)
    {
        vect3 &v=_operateurFerguson(x,j-offsetJ,m1);
        mat(offsetI+0,j)=v.X();
        mat(offsetI+1,j)=v.Y();
        mat(offsetI+2,j)=v.Z();
    }
}

void operateurDipolePotDer(const vect3 &r0, const vect3 &q,const mesh &inner_layer, vecteur &rhs, int offsetIdx)
{
    static analytiqueDipPotDer anaDPD;
    static integrateur<vect3> gauss;
    gauss.setOrdre(GaussOrder);

    for(int i=0;i<inner_layer.nbr_trg();i++)
    {
        anaDPD.init(inner_layer,i,q,r0);
        vect3 v=gauss.integre(anaDPD,inner_layer.trngl(i),inner_layer);
        rhs(inner_layer.trngl(i-offsetIdx).id1()+offsetIdx)+=v[0];
        rhs(inner_layer.trngl(i-offsetIdx).id2()+offsetIdx)+=v[1];
        rhs(inner_layer.trngl(i-offsetIdx).id3()+offsetIdx)+=v[2];
    }
}

void operateurDipolePot(const vect3 &r0, const vect3 &q, const mesh &inner_layer, vecteur &rhs, int offsetIdx)
{
    static analytiqueDipPot anaDP;
    static integrateur<double> gauss;
    gauss.setOrdre(GaussOrder);

    anaDP.init(q,r0);
    for(int i=offsetIdx;i<offsetIdx+inner_layer.nbr_trg();i++)
    {
        rhs(i)+=gauss.integre(anaDP,inner_layer.trngl(i-offsetIdx),inner_layer);
    }
}

// Grad wrt r0
void operateurDipolePotDerGrad(const vect3 &r0, const vect3 &q,const mesh &inner_layer, vecteur rhs[6], int offsetIdx)
{
    static analytiqueDipPotDerGrad anaDPD;
    static integrateur< vect3array<6> > gauss;
    gauss.setOrdre(GaussOrder);

    for(int i=0;i<inner_layer.nbr_trg();i++)
    {
        anaDPD.init(inner_layer,i,q,r0);
        vect3array<6> v=gauss.integre(anaDPD,inner_layer.trngl(i),inner_layer);
		for (int d=0;d<6;d++) { // derivatives wrt r0,q
			rhs[d](inner_layer.trngl(i-offsetIdx).id1()+offsetIdx)+=v(d)[0];
			rhs[d](inner_layer.trngl(i-offsetIdx).id2()+offsetIdx)+=v(d)[1];
			rhs[d](inner_layer.trngl(i-offsetIdx).id3()+offsetIdx)+=v(d)[2];
		}
    }
}

// Grad wrt r0,q
void operateurDipolePotGrad(const vect3 &r0, const vect3 &q, const mesh &inner_layer, vecteur rhs[6], int offsetIdx)
{
    static analytiqueDipPotGrad anaDP;
    static integrateur< vect3array<2> > gauss;
    gauss.setOrdre(GaussOrder);

    anaDP.init(q,r0);
    for(int i=offsetIdx;i<offsetIdx+inner_layer.nbr_trg();i++)
    {
        vect3array<2> v=gauss.integre(anaDP,inner_layer.trngl(i-offsetIdx),inner_layer);
		// grad_r0
		rhs[0](i)+=v(0).X();
		rhs[1](i)+=v(0).Y();
		rhs[2](i)+=v(0).Z();
		// grad_q
		rhs[3](i)+=v(1).X();
		rhs[4](i)+=v(1).Y();
		rhs[5](i)+=v(1).Z();
    }
}
