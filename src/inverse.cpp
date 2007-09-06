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
    string Option;
    matrice EegGainMatrix;
    sparse_matrice SmoothMatrix;
    fast_sparse_matrice fastSmoothMatrix;
    fast_sparse_matrice fastSmoothMatrix_t;
    vecteur AiVector;
    matrice RealEegData;
    matrice EstimatedSourcesData;
    double EegDataWeight;
    double SmoothWeight;
    string SmoothType;
    int MaxNbIter;
    double StoppingTol;
    matrice MegGainMatrix;
    double MegDataWeight;
    matrice RealMegData;
    
    // for use with EEG DATA
    if(!strcmp(argv[1],"-EEG"))
    {
        Option=string(argv[1]);
        EegGainMatrix.loadBin(argv[2]);
        SmoothMatrix.loadBin(argv[3]);
        fastSmoothMatrix=fast_sparse_matrice(SmoothMatrix);
        fastSmoothMatrix_t=fast_sparse_matrice(SmoothMatrix.transpose());
        AiVector.loadBin(argv[4]);
        RealEegData.loadTxt(argv[5]);
        EegDataWeight=atof(argv[7]);
        SmoothWeight=atof(argv[8]);
        SmoothType=string(argv[9]);
        MaxNbIter=atoi(argv[10]);
        StoppingTol=atof(argv[11]);
    }
    // for use with MEG DATA
    else if(!strcmp(argv[1],"-MEG"))
    {
        Option=string(argv[1]);
        MegGainMatrix.loadBin(argv[2]);
        SmoothMatrix.loadBin(argv[3]);
        fastSmoothMatrix=fast_sparse_matrice(SmoothMatrix);
        fastSmoothMatrix_t=fast_sparse_matrice(SmoothMatrix.transpose());
        AiVector.loadBin(argv[4]);
        RealMegData.loadTxt(argv[5]);
        MegDataWeight=atof(argv[7]);
        SmoothWeight=atof(argv[8]);
        SmoothType=string(argv[9]);
        MaxNbIter=atoi(argv[10]);
        StoppingTol=atof(argv[11]);
    }

    bool MEG,EEG;
    MEG=(Option==string("-MEG"));
    EEG=(Option==string("-EEG"));

    matrice &GainMatrix=MEG?MegGainMatrix:EegGainMatrix;
    matrice &data=MEG?RealMegData:RealEegData;
    double  DataWeight=MEG?MegDataWeight:EegDataWeight;
    size_t nT=data.ncol();
    EstimatedSourcesData=matrice(GainMatrix.ncol(),nT);
    vecteur v(EstimatedSourcesData.nlin());

    if(SmoothType==string("TV"))
    {
        for(size_t frame=0;frame<nT;frame++)
        {
            cout << ">> Frame " << frame+1 << endl;
            vecteur m_vec;
            m_vec.DangerousBuild(&data(0,frame),data.nlin());

            // ====================  initialization of source vector ===================== //
            if(frame==0) for(size_t i=0;i<v.size();i++) v(i)=0.0;//v(i)=1e-3*drandom(); // FIXME : add option for random init
            else for(size_t i=0;i<v.size();i++) v(i)=EstimatedSourcesData(i,frame-1);

            bool errorTest=true;
            double dtv=0.0;

            // ==========================  the inverse problem ========================== //
            int t;
            for(t=0;t<MaxNbIter && errorTest;t++)
            {
                vecteur gradtv=gentv(v,fastSmoothMatrix,fastSmoothMatrix_t,AiVector,&dtv);
                vecteur current_mes=GainMatrix*v;
                vecteur err_vec=current_mes-m_vec;
                vecteur graddata=GainMatrix.tmult(err_vec);
                vecteur Ggraddata=GainMatrix*graddata;

                double denom_data=Ggraddata*Ggraddata;
                double opt_step_data=-(Ggraddata*err_vec)/denom_data;
                vecteur grad=(-DataWeight)*graddata+(-SmoothWeight)*gradtv;
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

    bool Heat=SmoothType==string("HEAT");
    bool Mn=SmoothType==string("MN");

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
            m_vec.DangerousBuild(&data(0,frame),data.nlin());
            double alpha=SmoothWeight/DataWeight;

            //==========  initialization of source vector =======================//
            for(size_t i=0;i<v.size();i++) v(i)=0.0;//v(i)=1e-3*drandom(); // FIXME : add option for random init

            bool errorTest=true;
            double dtv=0;

            int t;
            for( t=0;t<MaxNbIter && errorTest;t++)
            {
                vecteur gradtv=gentv(v,fastSmoothMatrix,fastSmoothMatrix_t,AiVector,hess,&dtv,ftab[1],fptab[1],fpptab[1]);

                vecteur err_vec=GainMatrix*v-m_vec;
                vecteur graddata=GainMatrix.tmult(err_vec);
                vecteur grad=graddata+alpha*gradtv;

                LinOp &TIH=*(Heat?(LinOp*)new TvInverseHessian(GainMatrix,fastSmoothMatrix_t,hess,alpha):(LinOp*)new TikInverseHessian(GainMatrix,alpha));

                vecteur s(v.size()); s.set(0.0);
                for(size_t i=0;i<s.size();i++) s(i)=0.0;//s(i)=1e-3*drandom(); // FIXME : add option for random init

                MinRes2(TIH,-1.0*grad,s,MINRES_TOL);

                v=v+s;
                double move= (s*s)/(v*v);
                cout<<"Move = "<<move<<endl;
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
    // for use with EEG DATA
    if(!strcmp(argv[1],"-EEG"))
    {
        EstimatedSourcesData.saveTxt(argv[6]);
    }
    // for use with MEG DATA
    else if(!strcmp(argv[1],"-MEG"))
    {
        EstimatedSourcesData.saveTxt(argv[6]);
    }

    // Stop Chrono
    C.stop();
    C.dispEllapsed();

    return 0;
}

void getHelp(char** argv)
{
    cout << argv[0] <<" [-option] [filepaths...]" << endl << endl;

    cout << "-option :" << endl;
    cout << "   -EEG :   Compute the inverse for EEG " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            EegGainMatrix, SmoothMatrix, AiVector, RealEegData, EstimatedSourcesData, EegDataWeight, SmoothWeight, SmoothType, MaxNbIter, StoppingTol" << endl << endl;

    cout << "   -MEG :   Compute the inverse for MEG " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            MegGainMatrix, SmoothMatrix, AiVector, RealMegData, EstimatedSourcesData, MegDataWeight, SmoothWeight, SmoothType, MaxNbIter, StoppingTol" << endl << endl;

    exit(0);
}
