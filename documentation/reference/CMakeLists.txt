#
# Build the documentation
#

include(Documentation OPTIONAL)

if (BUILD_DOCUMENTATION)

    #
    # Configure the script and the doxyfile, then add target
    #
    # XXXX: This should not be done in this manner, or at least it should break when dot is not found.
    #       The way to do that is:
    #       find_package(Doxygen
    #                    REQUIRED dot
    #                    OPTIONAL_COMPONENTS gnuplot)
    #       if(NOT DOXYGEN_DOT_FOUND)
    #            message(WARNING "Graphviz doesn't seem to be installed.\
    #                             Doxygen will not be able to generate graphs.\
    #                             Consider installing this package.")
    #       endif(NOT DOXYGEN_DOT_FOUND)
    #
    if (NOT DOT_PATH)
        get_filename_component(DOT_PATH ${DOT} PATH)
    endif()


    if (NOT USE_CGAL)
        set(OpenMEEG_EXCLUDE_FROM_DOCUMENTATION
            "${OpenMEEG_EXCLUDE_FROM_DOCUMENTATION} ${CMAKE_SOURCE_DIR}/tools/cgal_mesh_create.cpp ${CMAKE_SOURCE_DIR}/tools/cgal_mesh_include.h")
    endif()
    configure_file(
        ${OpenMEEG_SOURCE_DIR}/documentation/reference/doxyfile.in
        ${OpenMEEG_BINARY_DIR}/documentation/reference/doxyfile)

    #   The ALL in the next directive is necessary to cope with a CMake bug !!!
    #   Otherwise the target is never activated.

    add_custom_target(OpenMEEG-documentation ALL)
    add_dependencies(OpenMEEG-documentation OpenMEEG-doxygen_outputs)
    add_custom_target(OpenMEEG-doxygen_outputs
                      COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile
                      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ DESTINATION share/doc/OpenMEEG/reference)

endif()
