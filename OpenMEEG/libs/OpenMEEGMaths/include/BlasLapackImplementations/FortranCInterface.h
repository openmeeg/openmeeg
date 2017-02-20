#pragma once

#if defined (WIN32)
    #define FC_GLOBAL(x,X) x ## _
#else
    #include <FC.h>
#endif
