#include "matrice.h"
#include "symmatrice.h"
#include "vecteur.h"
#include "sparse_matrice.h"
#include "fast_sparse_matrice.h"
#include "vect3.h"
#include "cpuChrono.h"

#include "om_utils.h"
#include "inverse.h"

using namespace std;

void getHelp(char** argv);

int main(int argc, char **argv)
{
    // Makefile for global usage
    if(argc==1)
    {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }

    if ((!strcmp(argv[1],"-h")) | (!strcmp(argv[1],"--help"))) getHelp(argv);

    // Start Chrono
    cpuChrono C;
    C.start();

    cout << endl << "| ------ " << argv[0] << " -------" << endl;
    for( int i = 1; i < argc; i += 1 )
    {
        cout << "| " << argv[i] << endl;
    }
    cout << "| -----------------------" << endl;

    // declaration of argument variables
    matrice GainMatrix;
    sparse_matrice SmoothMatrix;
    fast_sparse_matrice fastSmoothMatrix;
    fast_sparse_matrice fastSmoothMatrix_t;
    vecteur AiVector;
    matrice Data;
    matrice EstimatedSourcesData;
    double SmoothWeight;
    string SmoothType;
    int MaxNbIter;
    double StoppingTol;

    GainMatrix.loadBin(argv[1]);
    SmoothMatrix.loadBin(argv[2]);
    fastSmoothMatrix=fast_sparse_matrice(SmoothMatrix);
    fastSmoothMatrix_t=fast_sparse_matrice(SmoothMatrix.transpose());
    AiVector.loadBin(argv[3]);
    Data.loadTxt(argv[4]);
    SmoothWeight=atof(argv[6]);
    SmoothType=string(argv[7]);
    MaxNbIter=atoi(argv[8]);
    StoppingTol=atof(argv[9]);

    size_t nT=Data.ncol();
    EstimatedSourcesData=matrice(GainMatrix.ncol(),nT);
    vecteur v(EstimatedSourcesData.nlin());

    bool Heat = SmoothType==string("HEAT");
    bool Mn   = SmoothType==string("MN");
    bool Tv   = SmoothType==string("TV");

    if (!Tv && !Mn && !Heat) {
        std::cerr << "Unknown Smoothtype :  " << SmoothType << std::endl;
        std::cerr << "Should be HEAT , MN or TV" << std::endl;
        exit(1);
    }

    if(Tv)
    {
        for(size_t frame=0;frame<nT;frame++)
        {
            cout << ">> Frame " << frame+1 << endl;
            vecteur m_vec;
            m_vec.DangerousBuild(&Data(0,frame),Data.nlin());

            // ====================  initialization of source vector ===================== //
            if(frame==0) for(size_t i=0;i<v.size();i++) v(i)=0.0;//v(i)=1e-3*drandom(); // FIXME : add option for random init
            else for(size_t i=0;i<v.size();i++) v(i)=EstimatedSourcesData(i,frame-1);

            bool errorTest=true;
            double dtv=0.0;

            // ==========================  the inverse problem ========================== //
            int t;
            for(t=0;t<MaxNbIter && errorTest;t++)
            {
                vecteur gradtv = gentv(v,fastSmoothMatrix,fastSmoothMatrix_t,AiVector,&dtv);
                vecteur current_mes = GainMatrix*v;
                vecteur err_vec = current_mes-m_vec;
                vecteur graddata = GainMatrix.tmult(err_vec);
                vecteur Ggraddata = GainMatrix*graddata;

                double denom_data = Ggraddata*Ggraddata;
                double opt_step_data = -(Ggraddata*err_vec)/denom_data;
                vecteur grad = (-SmoothWeight)*gradtv - graddata;
                v=v+grad;
                double tol = sqrt((grad*grad)/(v*v));
                errorTest = tol>StoppingTol;
                if ((t%100)==0)
                    printf("TV= %f   Relative Error= %f   opt_step_data= %f   Tol= %f   Iter %d\n",dtv,(err_vec).norm()/m_vec.norm(),opt_step_data,tol,t);
            }
            //===========================================================================//
            for(size_t i=0;i<EstimatedSourcesData.nlin();i++) EstimatedSourcesData(i,frame)=v(i);
            m_vec.DangerousKill();
            std::cout << "Number of iterations = " << t << std::endl;
            cout << "Total Variation = " << dtv << endl;
        }
    }// TV SmoothType

    if(Heat || Mn)
    {
        // In both cases override maxnbiter number because the energy is quadratic and the minimum is found in one iteration
        // minor modifications should be done on the Newton algorithm for convex non-quadratic functions
        MaxNbIter=1;
        fast_sparse_matrice hess(fastSmoothMatrix);

        for(size_t frame=0;frame<nT;frame++)
        {
            cout << ">> Frame " << frame+1 << endl;
            vecteur m_vec;
            m_vec.DangerousBuild(&Data(0,frame),Data.nlin());

            //==========  initialization of source vector =======================//
            for(size_t i=0;i<v.size();i++) v(i)=0.0;//v(i)=1e-3*drandom(); // FIXME : add option for random init

            bool errorTest=true;
            double dtv=0;

            int t;
            for( t=0;t<MaxNbIter && errorTest;t++)
            {
                vecteur gradtv = gentv(v,fastSmoothMatrix,fastSmoothMatrix_t,AiVector,hess,&dtv,ftab[1],fptab[1],fpptab[1]);

                vecteur err_vec = GainMatrix*v-m_vec;
                vecteur graddata = GainMatrix.tmult(err_vec);
                vecteur grad = graddata+SmoothWeight*gradtv;

                LinOp &TIH = *( Heat ?
                                (LinOp*) new TvInverseHessian(GainMatrix,fastSmoothMatrix_t,hess,SmoothWeight) :
                                (LinOp*) new TikInverseHessian(GainMatrix,SmoothWeight) );

                vecteur s(v.size()); s.set(0.0);

                MinRes2(TIH,-1.0*grad,s,MINRES_TOL);

                v=v+s;
                double move = (s*s)/(v*v);
                std::cout << "Move = " << move << std::endl;
                if(move<StoppingTol) errorTest=false;

                if(Heat) delete (TvInverseHessian*)&TIH; else delete (TikInverseHessian*)&TIH;

            }// loop over t
            std::cout << "Number of iterations = " << t << std::endl;
            std::cout << "Relative Error = " << (GainMatrix*v-m_vec).norm()/m_vec.norm() << std::endl;
            for(size_t i=0;i<EstimatedSourcesData.nlin();i++) EstimatedSourcesData(i,frame)=v(i);
            m_vec.DangerousKill();

        }// loop over frame

    } // HEAT SmoothType

    // write output variables
    EstimatedSourcesData.saveTxt(argv[5]);

    // Stop Chrono
    C.stop();
    C.dispEllapsed();

    return 0;
}

void getHelp(char** argv)
{
    cout << argv[0] <<"[filepaths...]" << endl << endl;
    cout << "Compute the inverse for MEG/EEG " << endl;
    cout << "\tFilepaths are in order :" << endl;
    cout << "\tGainMatrix, SmoothMatrix, AiVector, RealData, EstimatedSourcesData, SmoothWeight, SmoothType, MaxNbIter, StoppingTol" << endl << endl;
    exit(0);
}
