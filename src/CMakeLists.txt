if (NOT HAVE_ISNORMAL_IN_NAMESPACE_STD)
    include(CheckSymbolExists)
    CHECK_SYMBOL_EXISTS(isnormal math.h HAVE_ISNORMAL_IN_MATH_H)
endif()

set(OPENMEEG_HEADERS
    analytics.h assemble.h cpuChrono.h danielsson.h DLLDefinesOpenMEEG.h domain.h forward.h gain.h geometry.h gmres.h integrator.h
    interface.h mesh.h om_utils.h operators.h options.h PropertiesSpecialized.h geometry_reader.h geometry_io.h sensors.h
    triangle.h Triangle_triangle_intersection.h vect3.h vertex.h 
#   These files are imported from another repository.
#   Please do not update them in this repository.
    DataTag.H FileExceptions.H GeometryExceptions.H Properties.H)

function(creta_library libname)
    add_library(${libname} ${LIB_TYPE} ${ARGN})
    set_target_properties(${libname} PROPERTIES
                             VERSION 1.1.0
                             SOVERSION 1
                             CLEAN_DIRECT_OUTPUT 1)
endfunction()

set(OpenMEEG_SOURCES 
    assembleFerguson.cpp assembleHeadMat.cpp assembleSourceMat.cpp assembleSensors.cpp domain.cpp triangle.cpp mesh.cpp interface.cpp
    danielsson.cpp geometry.cpp operators.cpp sensors.cpp cpuChrono.cpp ${OPENMEEG_HEADERS})

create_library(OpenMEEG ${OpenMEEG_SOURCES})
target_link_libraries(OpenMEEG OpenMEEGMaths ${OPENMEEG_LIBRARIES} ${LAPACK_LIBRARIES})
#create_library(OpenMEEGInteractive ${OpenMEEG_SOURCES})
#target_link_libraries(OpenMEEGInteractive OpenMEEGMathsInteractive ${OPENMEEG_LIBRARIES} ${LAPACK_LIBRARIES})

macro(OPENMEEG_EXECUTABLE EXECUTABLE)
    NEW_EXECUTABLE(${EXECUTABLE} ${ARGN} LIBRARIES OpenMEEG)
endmacro()

OPENMEEG_EXECUTABLE(om_minverser minverser.cpp)
OPENMEEG_EXECUTABLE(om_gain gain.cpp)
OPENMEEG_EXECUTABLE(om_forward forward.cpp)
OPENMEEG_EXECUTABLE(om_assemble assemble.cpp)

# ================
# = INSTALLATION =
# ================

install(TARGETS OpenMEEG om_assemble om_forward om_gain om_minverser
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

# install headers
install(FILES ${OPENMEEG_HEADERS}
        DESTINATION ${OPENMEEG_HEADER_INSTALLDIR} COMPONENT Development)

# ===========
# = TESTING =
# ===========

if (BUILD_TESTING)
    include(TestHead.txt)

    # those models should give same results !
    #TESTHEAD(0)
    #TESTHEAD(NNa0)
    TESTHEAD(1)
    TESTHEAD(NNa1)
    TESTHEAD(NNb1)
    TESTHEAD(NNc1)

    TESTHEAD(2)
    if (USE_VTK)
        TESTHEAD(NNa2)
        TESTHEAD(NNb2)
        TESTHEAD(NNc2)
    endif()

    option(TEST_HEAD3 "Run tests on Head 3" OFF)

    if (TEST_HEAD3)
        TESTHEAD(3)
        if (USE_VTK)
            TESTHEAD(NNa3)
            TESTHEAD(NNb3)
            TESTHEAD(NNc3)
        endif()
    endif()
endif()
