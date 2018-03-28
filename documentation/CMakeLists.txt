# # Require dot, treat the other components as optional
# find_package(Doxygen
#              REQUIRED dot) # XXX if dot module if not present, cmake should stop generatio
#              # OPTIONAL_COMPONENTS mscgen dia)

IF ( BUILD_DOCUMENTATION )
    SUBDIRS(reference)
ENDIF ( BUILD_DOCUMENTATION )

IF ( BUILD_TUTORIALS )
    SUBDIRS(tutorial)
ENDIF ( BUILD_TUTORIALS )
