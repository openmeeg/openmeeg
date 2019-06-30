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

#include <stdio.h>
#include <string>

#include <EEGsensors.h>
#include <MEGsensors.h>

using namespace OpenMEEG;

template <class S>
int test_sensors(const S& s) {

    size_t n = s.getNumberOfSensors();
    std::cout << "Number of sensors of S : " << n << std::endl;

    if (s.isEmpty())
        std::cout << "WARNING : empty sensors !" << std::endl;
    else {
        s.info();
        s.save("tmp.sensors");

        /**** test on copy constructor ****/
        S scopy("tmp.sensors");
        if (scopy.getNumberOfSensors() != n) {
            std::cout << "ERROR in copy from copy constructor : incorrect number of sensors" << std::endl;
            return -1;
        }

        scopy.info();
    }
    remove("tmp.sensors");
    return 0;
}

int main(const int argc, const char** argv) {
    if (argc < 2)
        return -1;

    std::string base_name = argv[1];
    EEGSensors eeg((base_name + ".eeg").c_str());
    MEGSensors meg((base_name + ".squids").c_str());
    return test_sensors(eeg) + test_sensors(meg);
}
