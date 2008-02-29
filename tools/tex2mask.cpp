/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
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

#include <matrice.h>
#include <symmatrice.h>
#include <vecteur.h>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

int main( int argc, char **argv)
{
    if ((argc != 4) && (argc != 5)) 
    {
        cout << "usage: "<< argv[0] <<" input_texture_file output_mask_txt_file min-thresh [max-thresh]" << endl;
        exit(1);
    }
    double min,max;
    bool has_max;
    if (argc == 4) 
    {
        min = atof(argv[3]);
    }
    if (argc == 5) 
    {
        min = atof(argv[3]);
        max = atof(argv[4]);
        has_max = true;
    }

    std::ifstream file(argv[1]);
    if (!file.is_open()) 
    {
        cout << "Problem reading : " << argv[1] << endl;
        return 1;
    }

    vector<string> tokens;

    std::string line;
    while (std::getline(file, line, '\n'))
    {
        if (!line.empty())
        {
            std::istringstream buffer(line);
            string tile;
            while (buffer >> tile)
            {
                tokens.push_back(tile);
            }
        }
    }

    vector<string>::iterator it = tokens.begin();

    if (*it != "ascii")
    {
        cout << "Not ascii format" << endl;
        return 1;
    }
    ++it;

    if (*it != "FLOAT")
    {
        cout << "Not FLOAT format but " << *it << endl;
        return 1;
    }
    ++it;

    int nb_time_steps = atoi((*it).c_str());
    cout << "Nb time steps : " << nb_time_steps << endl;
    ++it ;

    ++it ; // skipping header

    int nb_points = atoi((*it).c_str());
    cout << "Nb points : " << nb_points << endl;

    matrice m(nb_points,nb_time_steps);
    double val;

    for (int i = 0 ; i < nb_time_steps ; i++)
    {
        nb_points = atoi((*it).c_str());;
        cout << "Nb points at time " << i << " : " << nb_points << endl;
        ++it;
        for (int j = 0 ; j < nb_points ; j++)
        {
            val = (float) atof((*it).c_str());
            if(has_max)
            {
                if ((val >= min) && (val < max)) m(j,i) = 1;
                else m(j,i) = 0;
            }
            else {
                if (val >= min) m(j,i) = 1;
                else m(j,i) = 0;
            }
            ++it;
        }
        if (i != (nb_time_steps-1)) ++it;
    }

    m.saveTxt(argv[2]);

    return 0;
}
