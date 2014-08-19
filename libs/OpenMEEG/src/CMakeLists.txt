if (NOT HAVE_ISNORMAL_IN_NAMESPACE_STD)
    include(CheckSymbolExists)
    CHECK_SYMBOL_EXISTS(isnormal math.h HAVE_ISNORMAL_IN_MATH_H)
endif()

function(create_library libname)
    add_library(${libname} ${LIB_TYPE} ${ARGN})
    set_target_properties(${libname} PROPERTIES
                          VERSION 1.1.0
                          SOVERSION 1
                          CLEAN_DIRECT_OUTPUT 1)

    install(TARGETS ${libname}
            RUNTIME DESTINATION ${INSTALL_BIN_DIR}
            LIBRARY DESTINATION ${INSTALL_LIB_DIR}
            ARCHIVE DESTINATION ${INSTALL_LIB_DIR})
endfunction()

set(OpenMEEG_SOURCES 
    assembleFerguson.cpp assembleHeadMat.cpp assembleSourceMat.cpp assembleSensors.cpp domain.cpp mesh.cpp interface.cpp
    danielsson.cpp geometry.cpp operators.cpp sensors.cpp)

create_library(OpenMEEG ${OpenMEEG_SOURCES})
target_link_libraries(OpenMEEG PUBLIC OpenMEEGMaths PRIVATE ${OPENMEEG_LIBRARIES} ${matio_LIBRARIES} ${LAPACK_LIBRARIES})

#create_library(OpenMEEGInteractive ${OpenMEEG_SOURCES})
#target_link_libraries(OpenMEEGInteractive OpenMEEGMathsInteractive ${OPENMEEG_LIBRARIES} ${LAPACK_LIBRARIES})
