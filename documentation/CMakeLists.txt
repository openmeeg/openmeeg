# # Require dot, treat the other components as optional
# find_package(Doxygen
#              REQUIRED dot) # XXX if dot module if not present, cmake should stop generatio
#              # OPTIONAL_COMPONENTS mscgen dia)


if (BUILD_REFERENCE_DOC)# XXX: the right condition should be (BUILD_REFERENCE_DOC and DOXYGEN_FOUND) therefore, Doxygen should be found in documentation rather than its subdirectories
    subdirs(reference)
endif()
