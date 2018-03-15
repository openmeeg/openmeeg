#pragma once

#if defined (WIN32)
    #define FC_GLOBAL(x,X) x ## _
#elif defined (USE_VECLIB)
    // XXX : this is currently a hack as we did not have the mechanism to generate the FC.h
    #define FC_GLOBAL(name,NAME) NAME##_
#elif defined (USE_ATLAS)
    // XXX : this is currently a hack as we did not have the mechanism to generate the FC.h
    #define FC_GLOBAL(name,NAME) name##_
#else
    #include <FC.h>
#endif
