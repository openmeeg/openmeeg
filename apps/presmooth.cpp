/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

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

#include <symmatrix.h>
#include <vector.h>
#include <mesh.h>
#include <sparse_matrix.h>
#include <fast_sparse_matrix.h>

using namespace std;
using namespace OpenMEEG;

void getHelp(char** argv);

int main(int argc, char **argv) // TODO a quoi ça sert ?
{
    print_version(argv[0]);

    if (argc==1) {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }

    if ((!strcmp(argv[1],"-h")) | (!strcmp(argv[1],"--help"))) getHelp(argv);

    // declaration of argument variables
    Mesh SourceMesh;

    SourceMesh.load(argv[1]);

    SparseMatrix SmoothMatrix = SourceMesh.gradient();

    // write output variables
    SmoothMatrix.save(argv[2]);
    if (argc==4) {
        const Vector &AiVector = SourceMesh.areas();
        AiVector.save(argv[3]);
    }

    return 0;
}

void getHelp(char** argv)
{
    cout << argv[0] <<" [filepaths...]" << endl << endl;

    cout << "   Compute mesh gradient and triangles areas" << endl;
    cout << "   Filepaths are in order :" << endl;
    cout << "   SourceMesh, GradientMatrix (bin) [ , AreasVector (bin) ]" << endl << endl;

    exit(0);
}
