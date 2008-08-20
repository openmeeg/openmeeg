#ifndef _OpenMEEGMaths_DLLDEFINES_H_
#define _OpenMEEGMaths_DLLDEFINES_H_

/* Cmake will define OpenMEEGMaths_EXPORTS on Windows when it
configures to build a shared library. If you are going to use
another build system on windows or create the visual studio
projects by hand you need to define OpenMEEGMaths_EXPORTS when
building a DLL on windows.
*/

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
#endif

#if defined (WIN32)
  #if defined(OpenMEEGMaths_EXPORTS)
    #define  OPENMEEGMATHS_EXPORT __declspec(dllexport)
  #else
    #define  OPENMEEGMATHS_EXPORT __declspec(dllimport)
  #endif /* OpenMEEGMaths_EXPORTS */
#else /* defined (WIN32) */
 #define OPENMEEGMATHS_EXPORT
#endif

#endif /* _OpenMEEGMaths_DLLDEFINES_H_ */