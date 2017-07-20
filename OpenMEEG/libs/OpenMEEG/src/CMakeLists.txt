# =====================
# =    OpenMEEG Lib   =
# =====================

if (NOT HAVE_ISNORMAL_IN_NAMESPACE_STD)
    include(CheckSymbolExists)
    CHECK_SYMBOL_EXISTS(isnormal math.h HAVE_ISNORMAL_IN_MATH_H)
endif()

function(create_library libname)
    add_library(${libname} ${ARGN})
    if (LAPACK_DEFINITIONS)
        target_compile_definitions(${libname} PUBLIC ${LAPACK_DEFINITIONS})
    endif()
    set_target_properties(${libname} PROPERTIES
                          VERSION 1.1.0
                          SOVERSION 1
                          CLEAN_DIRECT_OUTPUT 1)

    generate_export_header(${libname}
        EXPORT_FILE_NAME ${libname}_Export.h
        STATIC_DEFINE ${libname}_BUILT_AS_STATIC)

    install(TARGETS ${libname}
            EXPORT  OpenMEEGTargets
            RUNTIME DESTINATION ${INSTALL_BIN_DIR}
            LIBRARY DESTINATION ${INSTALL_LIB_DIR}
            ARCHIVE DESTINATION ${INSTALL_LIB_DIR})

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${libname}_Export.h
            DESTINATION ${OpenMEEG_HEADER_INSTALLDIR} COMPONENT Development)
endfunction()

set(OpenMEEG_SOURCES 
    assembleFerguson.cpp assembleHeadMat.cpp assembleSourceMat.cpp assembleSensors.cpp domain.cpp
    mesh.cpp interface.cpp danielsson.cpp geometry.cpp operators.cpp sensors.cpp)

create_library(OpenMEEG ${OpenMEEG_SOURCES})
target_link_libraries(OpenMEEG PUBLIC OpenMEEGMaths ${LAPACK_LIBRARIES} ${VTK_LIBRARIES} ${GIFTI_LIBRARIES})
target_include_directories(OpenMEEG PUBLIC ${VTK_INCLUDE_DIRS} ${GIFTI_INCLUDE_DIRS})
