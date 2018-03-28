# # Require dot, treat the other components as optional
# find_package(Doxygen
#              REQUIRED dot) # XXX if dot module if not present, cmake should stop generatio
#              # OPTIONAL_COMPONENTS mscgen dia)


IF (BUILD_REFERENCE_DOC)# XXX: the right condition should be (BUILD_REFERENCE_DOC and DOXYGEN_FOUND) therefore, Doxygen should be found ind documentation rather than its subdirectories
    SUBDIRS(reference)
ENDIF (BUILD_REFERENCE_DOC)

IF (BUILD_TUTORIALS)
    SUBDIRS(tutorial)
ENDIF (BUILD_TUTORIALS)
