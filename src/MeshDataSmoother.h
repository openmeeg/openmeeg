#ifndef _H_MDS
#define _H_MDS

#include "mesh3.h"
#include <vector>
#include <set>
#include <algorithm>
#include "sparse_matrice.h"
#include "vecteur.h"

#define myFloat double

typedef sparse_matrice matType;

class MeshDataSmoother
{
private:
    std::vector< std::vector<int> >         SortedNeighbors;
    std::vector< myFloat >                  Parea;
    std::vector< std::vector<myFloat> >     cotanAlphaI;    // Pour l'isotrope
    std::vector< std::vector<myFloat> >     cotanBetaI;     // Pour l'isotrope
    std::vector< std::vector<myFloat> >     cotanThetaI;    // Pour l'anisotrope
    std::vector< std::vector<myFloat> >     cotanGammaI;    // Pour l'anisotrope
    std::vector< std::vector<myFloat> >     Tarea;
    std::vector< std::vector<Vect3> >       pip1p;          // Pour l'anisotrope
    std::vector< std::vector<Vect3> >       ppi;            // Pour l'anisotrope
    std::vector< myFloat >                  grad;           // Pour l'anisotrope
    myFloat*                                data;

    double dt;
    int n_iter;
    Mesh mesh;
    int nValues;

    void processMesh();
    void precalculate();
    bool aniso;

public:

    inline MeshDataSmoother(){}
    inline MeshDataSmoother( const char *filename, myFloat *dptr, int n_iter=0, const double &dt=0, bool aniso=false );
    inline ~MeshDataSmoother(){}

    inline void smooth(int n_time=0, const myFloat &dt=0, int* SubSourceIndexes=0, int SubSourceSize=0);
    inline vecteur getGradient(vecteur &x);
    inline myFloat computeMaxDt() const;
    inline Mesh& getMesh() {return mesh;}
    inline void setAniso ( const bool an ) {aniso = an;}

};

class MeshDataL1
{
private:
    Mesh mesh;
    int n;      // nombre de triangle
    int p;      // nombre de sommets
    matType *A; // matrice a calculer (3.n X p)
public:
    inline MeshDataL1();
    inline MeshDataL1(const char* filename);
    inline MeshDataL1(const Mesh &m);
    inline ~MeshDataL1();
    inline void computeMatrix();
    inline matType* getMatrix() {return A;}
};

class MeshDataL1Phi
{
private:
    Mesh mesh;
    int n;      // nombre de triangle^M
    int p;      // nombre de sommets^M
    matType *A; // matrice a calculer (3.n X p)^M
public:
    inline MeshDataL1Phi();
    inline MeshDataL1Phi(const char* filename);
    inline MeshDataL1Phi(const Mesh &m);
    inline ~MeshDataL1Phi();
    inline void computeMatrix(std::vector<double>* &pt);
    inline matType* getMatrix() {return A;}
};

class MeshDataL2
{
private:
    Mesh mesh;
    int n;      // nombre de sommets
    matType *A; //matrice a calculer
public:
    inline MeshDataL2();
    inline MeshDataL2(const char* filename);
    inline ~MeshDataL2();
    inline void computeMatrix();
    inline matType* getMatrix() {return A;}
};

inline vecteur chambolleSmooth( matType &Q, const vecteur &g, double lambda);

#define  EPSILON 1e-3
#include <fstream>

using namespace std;

inline int nextPointInCell ( const Mesh &mesh, int tid, int point)
{
    int pointIndexInCell=-1;
    for(int i=0;i<3;i++) if(mesh.getTrg(tid)[i]==point) pointIndexInCell = i;
    if(pointIndexInCell!=-1) return mesh.getTrg(tid)[(pointIndexInCell+1)%3];
    else return -1;
}

inline int prevPointInCell ( const Mesh &mesh, int tid, int point)
{
    int pointIndexInCell=-1;
    for(int i=0;i<3;i++) if(mesh.getTrg(tid)[i]==point) pointIndexInCell = i;
    if(pointIndexInCell!=-1) return mesh.getTrg(tid)[(pointIndexInCell+2)%3];
    else return -1;
}

inline MeshDataSmoother::MeshDataSmoother(const char *filename,myFloat *dptr_new,int n_iter_new, const double &dt_new, const bool an)
{
    mesh.load(filename);
    processMesh();
    this->data = dptr_new;
    this->dt = dt_new;
    this->n_iter = n_iter_new;
    nValues = mesh.nbPts();
    aniso = an;
}

inline void MeshDataSmoother::processMesh()
{
    intSet* pCells = mesh.getTrianglesForPoints(); // triangles to which a point belongs
    SortedNeighbors.resize(mesh.nbPts());
    std::vector<int> intersect(2);
    std::vector<int>::iterator intersect_it = intersect.begin();
    std::vector<int>::iterator intersect_it2;
    int pn,pnp1;
    grad = std::vector<myFloat>(mesh.nbPts());
    //ofstream fileout("c:\\outfile.txt");

    for(int i=0;i<mesh.nbPts();i++){

        // On remplit UnsortedNeighbors une premiere fois
        std::set<int> UnsortedNeighbors;
        for(std::set<int>::iterator it = pCells[i].begin(); it!=pCells[i].end(); it++){
            for(int ii=0;ii<3;ii++) UnsortedNeighbors.insert( mesh.getTrg(*it)[ii] );
        }
        SortedNeighbors[i]=std::vector<int>(UnsortedNeighbors.size());
        int cpt = 0;
        for(std::set<int>::iterator it=UnsortedNeighbors.begin(); it!=UnsortedNeighbors.end(); it++){
            if(*it!=i) SortedNeighbors[i].at(cpt++)=*it;
        }
        SortedNeighbors[i].resize(cpt);

        // On choisit
        pn = SortedNeighbors[i][0];

        // On prends le triangle de sorte que xPnPn+1 direct
        int n=(int)SortedNeighbors[i].size();
        cpt = 1;
        while(cpt<n){

            intersect_it2 = set_intersection(pCells[i].begin(), pCells[i].end(),
                pCells[pn].begin(), pCells[pn].end(),
                intersect_it);    // on a les deux triangles possibles
            int mydiff=(int)(intersect_it2-intersect_it);
            if(mydiff!=2) {
                cout<<"mydiff= "<<mydiff<<std::endl;
                for(intSet::iterator it = pCells[i].begin(); it!=pCells[i].end(); it++) cout<<(*it)<<"\t";
                assert(0);
            }

            if( prevPointInCell(mesh,intersect[0],pn)==i) pnp1 = nextPointInCell(mesh,intersect[0],pn);
            else pnp1 = nextPointInCell(mesh,intersect[1],pn);

            SortedNeighbors[i][cpt++]=pnp1;
            pn = pnp1;
        }
        SortedNeighbors[i].resize(n+1);
        SortedNeighbors[i][n]=SortedNeighbors[i][0]; // on fait un tableau circulaire.

    }    // loop over the points of the mesh
    // The sorted vincinities are now computed

    Tarea.resize(mesh.nbPts());
    Parea.resize(mesh.nbPts());
    cotanGammaI.resize(mesh.nbPts());
    cotanThetaI.resize(mesh.nbPts());
    cotanAlphaI.resize(mesh.nbPts());
    cotanBetaI.resize(mesh.nbPts());
    pip1p.resize(mesh.nbPts());
    ppi.resize(mesh.nbPts());


    for(int j=0;j<mesh.nbPts();j++)
    {
        int n=(int)SortedNeighbors[j].size()-1;
        Tarea[j]=std::vector<myFloat>(n);
        cotanGammaI[j]=std::vector<myFloat>(n);
        cotanThetaI[j]=std::vector<myFloat>(n);
        cotanAlphaI[j]=std::vector<myFloat>(n);
        cotanBetaI[j]=std::vector<myFloat>(n);
        pip1p[j]=std::vector<Vect3>(n);
        ppi[j]=std::vector<Vect3>(n);

        Parea[j]=0;
        for( int i=0; i<n; i++ )
        {
            const Vect3 &p = mesh.getPt(j);
            const Vect3 &pi = mesh.getPt(SortedNeighbors[j][i]);
            const Vect3 &pip1 = mesh.getPt(SortedNeighbors[j][i+1]);
            Tarea[j][i]=1.0 / max(( ( pip1 - p ) ^ ( pi - p )).norme(),EPSILON);
            Parea[j]+=0.5/Tarea[j][i];
            cotanGammaI[j][i]= Tarea[j][i]*( ( pi - pip1 ) * ( pi - p ) );
            cotanThetaI[j][i]= Tarea[j][i]*( ( pip1 - p ) * ( pip1 - pi ) );
            pip1p[j][i]=p-pip1;
            ppi[j][i]=pi-p;
            cotanAlphaI[j][i]= Tarea[j][i]*( ( pi - p ) * ( pi - pip1 ) );
            cotanBetaI[j][i]= Tarea[j][i]*( ( pip1 - p ) * ( pi - pip1 ) );
            }
        Parea[j]=1.0/max(Parea[j],EPSILON);
    }

    delete[] pCells;

}

inline void MeshDataSmoother::smooth(int N_iter, const myFloat &Dt, int* SubSourceIndexes, int SubSourceSize)
{
    if(N_iter!=0) this->n_iter = N_iter;
    if(Dt!=0)      this->dt = Dt;

    if(!SubSourceSize)
    {
        if(aniso)
        {
            for(int iter=0;iter<n_iter;iter++)
            {

                for(int iv=0;iv<nValues;iv++)
                {
                    myFloat myArea;
                    myArea = 0;
                    grad[iv] = 0;

                    for(int i=0;i<(int)SortedNeighbors[iv].size()-1;i++)
                    {
                        myFloat &Fp = data[iv];
                        myFloat &Fpi = data[SortedNeighbors[iv][i]];
                        myFloat &Fpip1 = data[SortedNeighbors[iv][i+1]];

                        myFloat IntrGradNorm = Tarea[iv][i]* (myFloat)( (Fpi-Fp)*(pip1p[iv][i]) + (Fpip1-Fp)*(ppi[iv][i])).norme();

                        grad[iv]+= 1/max(IntrGradNorm,EPSILON)* ( cotanThetaI[iv][i]*(Fpi-Fp)+cotanGammaI[iv][i]*(Fpip1-Fp));
                    }

                    grad[iv]*=Parea[iv];

                }

                for(int iv=0;iv<nValues;iv++) {data[iv]+=dt*grad[iv];
                }

            }
        }
        else{

            for(int iter=0;iter<n_iter;iter++)
            {

                for(int iv=0;iv<nValues;iv++){
                    myFloat myArea;
                    myArea = 0;
                    grad[iv] = 0;

                    for(int i = 0;i<(int)SortedNeighbors[iv].size()-1;i++)
                    {
                        myFloat &Fp = data[iv];
                        myFloat &Fpi = data[SortedNeighbors[iv][i]];
                        myFloat &Fpip1 = data[SortedNeighbors[iv][i+1]];

                        grad[iv]+= -cotanGammaI[iv][i]*(Fp-Fpip1) + cotanThetaI[iv][i]*(Fpi-Fp);
                    }

                    grad[iv]*=Parea[iv];

                }

                for(int iv=0;iv<nValues;iv++) {
                    data[iv]+=dt*grad[iv];
                }
            }

        }

    }
    else // if subsourcesize!=0
    {
        if(aniso)
        {
            for(int iter=0;iter<n_iter;iter++)
            {

                for(int iv=0;iv<SubSourceSize;iv++)
                {
                    int iiv = SubSourceIndexes[iv];
                    myFloat myArea;
                    myArea = 0;
                    grad[iiv] = 0;

                    for(int i=0;i<(int)SortedNeighbors[iiv].size()-1;i++)
                    {
                        myFloat &Fp = data[iiv];
                        myFloat &Fpi = data[SortedNeighbors[iiv][i]];
                        myFloat &Fpip1 = data[SortedNeighbors[iiv][i+1]];

                        myFloat IntrGradNorm = Tarea[iiv][i]* (myFloat)( (Fpi-Fp)*(pip1p[iiv][i]) + (Fpip1-Fp)*(ppi[iiv][i])).norme();

                        grad[iiv]+= 1/max(IntrGradNorm,EPSILON)* ( cotanThetaI[iiv][i]*(Fpi-Fp)+cotanGammaI[iiv][i]*(Fpip1-Fp));
                    }

                    grad[iiv]*=Parea[iiv];

                }

                for(int iv=0;iv<SubSourceSize;iv++) {data[SubSourceIndexes[iv]]+=dt*grad[SubSourceIndexes[iv]]; //moy+=fabs(dt*grad[iv]);
                }

            }
        }
        else
        {
            //double moy=0;

            for(int iter=0;iter<n_iter;iter++)
            {

                for(int iv=0;iv<SubSourceSize;iv++)
                {
                    int iiv = SubSourceIndexes[iv];
                    myFloat myArea;
                    myArea = 0;
                    grad[iiv] = 0;

                    for(int i=0;i<(int)SortedNeighbors[iiv].size()-1;i++)
                    {
                        myFloat &Fp = data[iiv];
                        myFloat &Fpi = data[SortedNeighbors[iiv][i]];
                        myFloat &Fpip1 = data[SortedNeighbors[iiv][i+1]];

                        grad[iiv] += -cotanGammaI[iiv][i]*(Fp-Fpip1) + cotanThetaI[iiv][i]*(Fpi-Fp);
                    }

                    grad[iiv]*=Parea[iiv];

                }

                for(int iv=0;iv<SubSourceSize;iv++)
                {
                    data[SubSourceIndexes[iv]]+=dt*grad[SubSourceIndexes[iv]];
                }

            }

        }


    }
}

inline vecteur MeshDataSmoother::getGradient(vecteur &x)
{
    myFloat *oldData = this->data;
    this->data=(myFloat*)& x(0);
    vecteur result(x.size());

    if(aniso)
    {
        for(int iv=0;iv<nValues;iv++)
        {
            myFloat myArea;
            myArea = 0;
            grad[iv] = 0;

            for(int i = 0;i<(int)SortedNeighbors[iv].size()-1;i++)
            {
                myFloat &Fp = data[iv];
                myFloat &Fpi = data[SortedNeighbors[iv][i]];
                myFloat &Fpip1 = data[SortedNeighbors[iv][i+1]];

                myFloat IntrGradNorm = Tarea[iv][i]* (myFloat)( (Fpi-Fp)*(pip1p[iv][i]) + (Fpip1-Fp)*(ppi[iv][i])).norme();

                grad[iv]+= 1/max(IntrGradNorm,EPSILON)* ( cotanThetaI[iv][i]*(Fpi-Fp)+cotanGammaI[iv][i]*(Fpip1-Fp));
            }

            grad[iv]*=Parea[iv];

        }

        for(int iv=0;iv<nValues;iv++)
        {
            result(iv)=-grad[iv]; //moy+=fabs(dt*grad[iv]);
        }

    }
    else{

        for(int iv=0;iv<nValues;iv++)
        {
            myFloat myArea;
            myArea = 0;
            grad[iv] = 0;

            for(int i = 0;i<(int)SortedNeighbors[iv].size()-1;i++)
            {
                myFloat &Fp = data[iv];
                myFloat &Fpi = data[SortedNeighbors[iv][i]];
                myFloat &Fpip1 = data[SortedNeighbors[iv][i+1]];

                grad[iv]+= -cotanGammaI[iv][i]*(Fp-Fpip1) + cotanThetaI[iv][i]*(Fpi-Fp);
            }

            grad[iv]*=Parea[iv];

        }

        for(int iv=0;iv<nValues;iv++) {result(iv)=-grad[iv];// moy+=fabs(dt*grad[iv]);
        }

    }

    this->data = oldData;

    return result;
}

inline myFloat MeshDataSmoother::computeMaxDt() const
{
    myFloat maxDt = 1e99;

    for(int iv=0;iv<nValues;iv++)
    {
        myFloat myArea;
        myArea = 0;
        myFloat gradient = 0;
        myFloat minFpi = 1e99;
        myFloat maxFpi = -1e99;

        myFloat &Fp = data[iv];

        for(int i=0;i<(int)SortedNeighbors[iv].size()-1;i++)
        {
            myFloat &Fpi = data[SortedNeighbors[iv][i]];
            myFloat &Fpip1 = data[SortedNeighbors[iv][i+1]];

            minFpi = min(Fpi,minFpi);
            maxFpi = max(Fpi,maxFpi);

            gradient+= cotanGammaI[iv][i]*(Fp-Fpip1) - cotanThetaI[iv][i]*(Fpi-Fp);
        }

        minFpi = abs(Fp-minFpi);
        maxFpi = abs(Fp-maxFpi);

        gradient*=Parea[iv];

        if(gradient!=0) maxDt = min(maxDt, min(minFpi,maxFpi)/abs(gradient));

        assert(maxDt!=0);

    }
    return maxDt;

}

inline Vect3 P1Vector( const Vect3 &p0, const Vect3 &p1, const Vect3 &p2, const int idx )
{
    assert(idx>-1 && idx<3);
    int i = idx+1;
    Vect3 pts[5] = {p2,p0,p1,p2,p0};
    Vect3 ret(0,0,0);
    Vect3 pim1pi = pts[i]-pts[i-1];
    Vect3 pim1pip1 = pts[i+1]-pts[i-1];
    Vect3 pim1H = ( (1.0/pim1pip1.norme2()) * ( pim1pi*pim1pip1 ) ) *pim1pip1;
    Vect3 piH = pim1H-pim1pi;
    ret = -1.0/piH.norme2()*piH;

    return ret;
}

inline double aire( const Vect3 &p0, const Vect3 &p1, const Vect3 &p2)
{
    Vect3 p0p1 = p1-p0;
    Vect3 p0p2 = p2-p0;

    return 0.5*(p0p1^p0p2).norme();
}


inline MeshDataL1::MeshDataL1(const char* filename)
{
    mesh.load(filename);
    n = mesh.nbTrgs();
    p = mesh.nbPts();
    A = new matType(3*n,p);
}

inline MeshDataL1::MeshDataL1(const Mesh &m)
{
    mesh = m;
    n = mesh.nbTrgs();
    p = mesh.nbPts();
    A = new matType(3*n,p);
}

inline MeshDataL1::~MeshDataL1()
{
    delete A;
}

inline void MeshDataL1::computeMatrix()
{
    //parcours des triangles
    for(int t=0;t<mesh.nbTrgs();t++)
    {
        const Triangle& trg = mesh.getTrg(t);
        Vect3 pts[3] = {mesh.getPt(trg[0]) ,mesh.getPt(trg[1]), mesh.getPt(trg[2])};
        Vect3 grads[3];
        double K = aire(pts[0], pts[1], pts[2]);
        for(int i=0;i<3;i++) grads[i]=P1Vector(pts[0], pts[1], pts[2], i);

        for(int i=0;i<3;i++) {
            for(int j=0;j<3;j++)
            {
                (*A)(3*t+i,trg[j]) = K*grads[j][i];
            }
        }
    }
}


inline MeshDataL1Phi::MeshDataL1Phi(const char* filename)
{
    mesh.load(filename);
    n = mesh.nbTrgs();
    p = mesh.nbPts();
    A = new matType(3*n,p);
}

inline MeshDataL1Phi::MeshDataL1Phi(const Mesh &m)
{
    mesh = m;
    n = mesh.nbTrgs();
    p = mesh.nbPts();
    A = new matType(3*n,p);
}

inline MeshDataL1Phi::~MeshDataL1Phi()
{
    delete A;
}

inline void MeshDataL1Phi::computeMatrix(std::vector<double>* &pt)
{
    pt = new std::vector<double>;
    // loop on triangles
    for(int t=0;t<mesh.nbTrgs();t++)
    {
        const Triangle& trg = mesh.getTrg(t);
        Vect3 pts[3] = {mesh.getPt(trg[0]) ,mesh.getPt(trg[1]), mesh.getPt(trg[2])};
        Vect3 grads[3];
        pt->push_back(aire(pts[0], pts[1], pts[2]));
        for(int i=0;i<3;i++) grads[i]=P1Vector(pts[0], pts[1], pts[2], i);

        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                (*A)(3*t+i,trg[j]) = grads[j][i];
            }
        }
    }
}

inline MeshDataL2::MeshDataL2(const char* filename)
{
    mesh.load(filename);
    n = mesh.nbPts();
    A = new matType(n,n);
}

inline MeshDataL2::~MeshDataL2()
{
    delete A;
}

inline void MeshDataL2::computeMatrix()
{
    // loop on triangles
    for(int t=0;t<mesh.nbTrgs();t++)
    {
        const Triangle& trg = mesh.getTrg(t);
        Vect3 pts[3] = {mesh.getPt(trg[0]) ,mesh.getPt(trg[1]), mesh.getPt(trg[2])};
        Vect3 grads[3];
        double K = aire(pts[0], pts[1], pts[2]);
        for(int i=0;i<3;i++) grads[i]=P1Vector(pts[0], pts[1], pts[2], i);

        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                (*A)(trg[i],trg[j])+=K*(grads[i]*grads[j]); // suppose symetry in matrix addresses
            }
        }
    }

    std::cout<<"MeshDataL2::computeMatrix()  WARNING"<<std::endl;
}

inline double norm3 ( const vecteur& v, size_t i)
{
    return sqrt(v(3*i)*v(3*i)+v(3*i+1)*v(3*i+1)+v(3*i+2)*v(3*i+2));
}

inline vecteur chambolleSmooth( matType &Q, const vecteur &g, double lambda)
{
    static const double tau = 0.01;
    static const double tol = 1e-3;

    vecteur p(Q.nlin()); p.set(0);
    vecteur p_old(Q.nlin()); p.set(0);
    vecteur q(Q.nlin());
    matType Qstar = Q.transpose();

    double delta;

    do
    {
        p_old = p.duplicate();
        vecteur qtemp = Qstar*p-g*(1.0/lambda);
        q = Q*qtemp;

        for(size_t i=0;i<Q.nlin()/3;i++)
        {
            double coef = 1.0/(1+tau*norm3(q,i));
            double coefXtau = tau*coef;
            for(int k=0;k<3;k++)
            {
                size_t h = 3*i+k;
                p(h)*=coef;
                p(h)-=coefXtau*q(h);
            }
        }
        delta = (p_old-p).norm()/p_old.norm();
        //cout<<"delta="<<delta<<std::endl;
    }
    while(delta>tol);

    return (Qstar*p)*lambda;
}


#endif

