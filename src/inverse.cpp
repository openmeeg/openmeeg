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

#include "cpuChrono.h"
#include "om_utils.h"
#include "inversers.h"

using namespace std;
using namespace OpenMEEG;

void getHelp(char** argv);

int main(int argc, char **argv)
{
    print_version(argv[0]);

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
    
    disp_argv(argc,argv);
    
    // declaration of argument variables
    Matrix GainMatrix;
    SparseMatrix SmoothMatrix;
    Vector AiVector;
    Matrix Data;
    Matrix EstimatedSourcesData;
    double SmoothWeight;
    string SmoothType;
    
    GainMatrix.loadBin(argv[1]);
    SmoothMatrix.loadBin(argv[2]);
    AiVector.loadBin(argv[3]);
    Data.loadTxt(argv[4]);
    SmoothWeight = atof(argv[6]);
    SmoothType   = string(argv[7]);
    
    bool Heat = SmoothType==string("HEAT");
    bool Mn   = SmoothType==string("MN");
    bool IMn  = SmoothType==string("IMN");
    bool WMn  = SmoothType==string("WMN");
    bool Tv   = SmoothType==string("TV");
    
    if (!Tv && !Mn && !IMn && !Heat && !WMn) {
        std::cerr << "Unknown Smoothtype :  " << SmoothType << std::endl;
        std::cerr << "Should be HEAT, IMN, MN or TV" << std::endl;
        exit(1);
    }
    
    if(Tv)
    {
        size_t MaxNbIter   = (size_t) atoi(argv[8]);
        double StoppingTol = atof(argv[9]);
    
        TV_inverse EstimatedSourcesData(Data,GainMatrix,SmoothMatrix,AiVector,SmoothWeight,MaxNbIter,StoppingTol);
        EstimatedSourcesData.saveTxt(argv[5]);
    }
    
    if(Mn)
    {
        MN_inverse EstimatedSourcesData(Data,GainMatrix,SmoothWeight);
        EstimatedSourcesData.saveTxt(argv[5]);
    }
    
    if(IMn)
    {
        IMN_inverse EstimatedSourcesData(Data,GainMatrix,SmoothWeight);
        EstimatedSourcesData.saveTxt(argv[5]);
    }
    
    if(WMn)
    {
        WMN_inverse EstimatedSourcesData(Data,GainMatrix,SmoothWeight);
        EstimatedSourcesData.saveTxt(argv[5]);
    }
    
    if(Heat)
    {
        HEAT_inverse EstimatedSourcesData(Data,GainMatrix,SmoothMatrix,SmoothWeight);
        EstimatedSourcesData.saveTxt(argv[5]);
    }
    
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
