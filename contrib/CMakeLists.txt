IF (USE_MATIO)

    #   Duplicate zlib discovery here in order not to patch Matio...

    OPTION(USE_SYSTEM_ZLIB "Use installed zlib if possible" ON)
    MARK_AS_ADVANCED(USE_SYSTEM_ZLIB)

    IF (USE_SYSTEM_ZLIB)

        #   Always impose a dynamic libz because it is linked (via matio) to Matlib which needs to
        #   be dynamic. Unfortunately, old linux systems offered libz.a which were not position independant.
        #   This is the modification of matio for OpenMEEG.

        IF (UNIX AND NOT APPLE)
            FIND_LIBRARY(SYSTEM_ZLIB z NAMES libz.so)
        ELSE()
            FIND_LIBRARY(SYSTEM_ZLIB z)
        ENDIF()
    ENDIF()

    SUBDIRS(matio)
ENDIF()

IF (NOT USE_MKL AND NEED_CLAPACK)
    SUBDIRS(clapack-3.2.1-CMAKE)
ENDIF()
