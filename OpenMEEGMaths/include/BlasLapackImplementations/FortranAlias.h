#pragma once

// Accelerate exports these as e.g. "_dgemm$NEWLAPACK", unspellable as a C
// identifier; without the asm label we would silently bind to Apple's frozen,
// deprecated LAPACK 3.2.1 symbols instead.
#if defined(USE_ACCELERATE)
    #define FC_ALIAS(x) __asm("_" #x "$NEWLAPACK")
#else
    #define FC_ALIAS(x)
#endif
