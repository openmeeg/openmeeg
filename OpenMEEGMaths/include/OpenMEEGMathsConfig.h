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

#pragma once

//  cmake configuration.

#include <OpenMEEGConfigure.h>
#include <OpenMEEGMaths_Export.h>

// specially for windows
#if WIN32
    #pragma inline_recursion (on)
    #pragma inline_depth (255) // MSVC static build with MKL cause LNK2019 error
    #pragma warning( disable : 4530)    //MSVC standard library can't be inlined
    #pragma warning( disable : 4996)    //MSVC warning C4996: declared deprecated
    #if defined(_MSC_VER)
        // Enable MSVC compiler warning messages that are useful but off by default.
        # pragma warning ( default : 4263 ) /* no override, call convention differs */
        // Disable MSVC compiler warning messages that often occur in valid code.
        # pragma warning ( disable : 4097 ) /* typedef is synonym for class */
        # pragma warning ( disable : 4127 ) /* conditional expression is constant */
        # pragma warning ( disable : 4244 ) /* possible loss in conversion */
        # pragma warning ( disable : 4251 ) /* missing DLL-interface */
        # pragma warning ( disable : 4305 ) /* truncation from type1 to type2 */
        # pragma warning ( disable : 4309 ) /* truncation of constant value */
        # pragma warning ( disable : 4514 ) /* unreferenced inline function */
        # pragma warning ( disable : 4706 ) /* assignment in conditional expression */
        # pragma warning ( disable : 4710 ) /* function not inlined */
        # pragma warning ( disable : 4786 ) /* identifier truncated in debug info */
        # pragma warning ( disable : 4244 ) /* possible loss of data ('float' to 'mat_uint32_t') */
        # pragma warning ( disable : 4267 ) /* possible loss of data (size_t to int) */
    #endif
#endif

//  Blas/Lapack configuration

#if defined(MKL_ILP64) && defined(USE_MKL)
typedef long long int BLAS_INT;
#else
typedef int BLAS_INT;
#endif

#if defined(USE_LAPACK)
#include <BlasLapackImplementations/OpenMEEGMathsBlasLapackConfig.h>
#elif defined(USE_MKL)
#include <BlasLapackImplementations/OpenMEEGMathsMKLConfig.h>
#elif defined(USE_ATLAS)
#include <BlasLapackImplementations/OpenMEEGMathsAtlasConfig.h>
#elif defined(USE_OPENBLAS)
#include <BlasLapackImplementations/OpenMEEGMathsOpenBLASConfig.h>
#elif defined(USE_VECLIB)
#include <BlasLapackImplementations/OpenMEEGMathsvecLibConfig.h>
#else
#warning "No blas/lapack implementation selected."
#endif

//#define inline __forceinline
//#define inline __attribute__((always_inline))
//#define inline __attribute__((weak)) inline


#define DPOTF2 LAPACK(dpotf2,DPOTF2)
#define DSPEVD LAPACK(dspevd,DSPEVD)
