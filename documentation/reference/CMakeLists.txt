#
# Build the documentation
#
INCLUDE (Documentation OPTIONAL)

IF (BUILD_DOCUMENTATION)

    #
    # Configure the script and the doxyfile, then add target
    #
    IF(NOT DOT_PATH)
        GET_FILENAME_COMPONENT(DOT_PATH ${DOT} PATH)
    ENDIF(NOT DOT_PATH)

    CONFIGURE_FILE(
        ${OpenMEEG_SOURCE_DIR}/documentation/reference/doxyfile.in
        ${OpenMEEG_BINARY_DIR}/documentation/reference/doxyfile)

    ADD_CUSTOM_COMMAND(OUTPUT html/index.html
                       COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile
                       DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)

    INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION share/doc/OpenMEEG/reference)

ENDIF (BUILD_DOCUMENTATION)
