INCLUDE (Documentation OPTIONAL)

ADD_CUSTOM_TARGET(documentation ALL)

ADD_DEPENDENCIES(documentation doxygen_outputs)

ADD_CUSTOM_TARGET(doxygen_outputs
                  COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config)

INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ DESTINATION share/doc/OpenMEEG/reference)
