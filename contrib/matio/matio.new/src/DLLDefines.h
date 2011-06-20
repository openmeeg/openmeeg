// Contents of DLLDefines.h
#ifndef MATIO_DLLDEFINES_H_
#define MATIO_DLLDEFINES_H_

/* Cmake will define MATIO_EXPORTS on Windows when it
configures to build a shared library. If you are going to use
another build system on windows or create the visual studio
projects by hand you need to define MyLibrary_EXPORTS when
building a DLL on windows.
*/
// We are using the Visual Studio Compiler and building Shared libraries

#if defined (WIN32) 
  #if defined(MATIO_INTERNAL)
    #define  MATIO_EXPORT __declspec(dllexport)
  #else
    #define  MATIO_EXPORT __declspec(dllimport)
  #endif /* matio_EXPORTS */
#else /* defined (WIN32) */
 #define MATIO_EXPORT
#endif

#endif /* MATIO_DLLDEFINES_H_ */
