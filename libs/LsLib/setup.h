/*! \addtogroup reelGroup reel and associated functions
    \brief Here is a list of functions close to reel
*/

/*! \file
    \brief file defining the precision of the type \c reel

    \c reel is either \c float or \c double
*/
#ifndef _SETUP_H
#define _SETUP_H

/*! \ingroup reelGroup */
/*!@{*/

/*!    \brief  flag to set the type of \c reel and the value of \c REEL_MAX
        - 0:
            - \c reel = \c double
            - \c REEL_MAX = \c DBL_MAX
        - 1:
            - \c reel = \c float
            - \c REEL_MAX = \c FLT_MAX


*/
#define USE_FLOATS 0

#include <float.h>
#include <limits.h>

#if USE_FLOATS
typedef float reel;         //!< definition of reel (\c USE_FLOATS = 1)
#define REEL_MAX FLT_MAX    //!< definition of REEL_MAX (\ USE_FLOATS = 1)
#else
typedef double reel;        //!< definition of reel (USE_FLOATS = 0)
#define REEL_MAX DBL_MAX    //!< definition of REEL_MAX (USE_FLOATS = 0)
#endif

/*!@}*/
#endif
