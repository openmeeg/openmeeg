#
# Build the documentation
#
INCLUDE (${CMAKE_ROOT}/Modules/Documentation.cmake OPTIONAL)

IF (BUILD_DOCUMENTATION)

    #
    # Configure the script and the doxyfile, then add target
    #
    IF(NOT DOT_PATH)
        GET_FILENAME_COMPONENT(DOT_PATH ${DOT} PATH)
    ENDIF(NOT DOT_PATH)

    CONFIGURE_FILE(
        ${OpenMEEG_SOURCE_DIR}/documentation/doxyfile.in
        ${OpenMEEG_BINARY_DIR}/documentation/doxyfile)

    ADD_CUSTOM_TARGET(documentation
        ${DOXYGEN_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/documentation/doxyfile)

ENDIF (BUILD_DOCUMENTATION)
