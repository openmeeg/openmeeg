// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

//  cmake configuration.

#include <OpenMEEGConfigure.h>
#include <OpenMEEGMaths_Export.h>

// specially for windows
#ifdef WIN32
    #pragma inline_recursion(on)
    #pragma inline_depth(255) // MSVC static build with MKL cause LNK2019 error
    #if defined(_MSC_VER)
        // Enable MSVC compiler warning messages that are useful but off by default.
        #pragma warning(default : 4263) /* no override, call convention differs */
        // Disable MSVC compiler warning messages that often occur in valid code.
        #pragma warning(disable : 4097) /* typedef is synonym for class */
        #pragma warning(disable : 4127) /* conditional expression is constant */
        #pragma warning(disable : 4244) /* possible loss in conversion */
        #pragma warning(disable : 4251) /* missing DLL-interface */
        #pragma warning(disable : 4305) /* truncation from type1 to type2 */
        #pragma warning(disable : 4309) /* truncation of constant value */
        #pragma warning(disable : 4514) /* unreferenced inline function */
        #pragma warning(disable : 4706) /* assignment in conditional expression */
        #pragma warning(disable : 4710) /* function not inlined */
        #pragma warning(disable : 4786) /* identifier truncated in debug info */
        #pragma warning(disable : 4244) /* possible loss of data ('float' to 'mat_uint32_t') */
        #pragma warning(disable : 4267) /* possible loss of data (size_t to int) */
    #endif
#endif

//  Blas/Lapack configuration

#if defined(USE_LAPACK)
#include <BlasLapackImplementations/OpenMEEGMathsBlasLapackConfig.h>
#elif defined(USE_MKL)
#include <BlasLapackImplementations/OpenMEEGMathsMKLConfig.h>
#elif defined(USE_ATLAS)
#include <BlasLapackImplementations/OpenMEEGMathsAtlasConfig.h>
#elif defined(USE_OPENBLAS)
#include <BlasLapackImplementations/OpenMEEGMathsOpenBLASConfig.h>
#elif defined(USE_FLEXIBLAS)
#include <BlasLapackImplementations/OpenMEEGMathsFlexiBLASConfig.h>
#else
#warning "No blas/lapack implementation selected."
#endif

#define DPOTF2 LAPACK(dpotf2,DPOTF2)
#define DSPEVD LAPACK(dspevd,DSPEVD)
