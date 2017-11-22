macro(OPENMEEG_EXECUTABLE EXECUTABLE)
    # NEW_EXECUTABLE(${EXECUTABLE} ${ARGN} LIBRARIES OpenMEEG ${LAPACK_LIBRARIES} ${VTK_LIBRARIES} ${GIFTI_LIBRARIES})
    NEW_EXECUTABLE(${EXECUTABLE} ${ARGN} LIBRARIES OpenMEEG ${VTK_LIBRARIES} ${GIFTI_LIBRARIES})
endmacro()

OPENMEEG_EXECUTABLE(om_minverser minverser.cpp)
OPENMEEG_EXECUTABLE(om_gain gain.cpp)
OPENMEEG_EXECUTABLE(om_forward forward.cpp)
OPENMEEG_EXECUTABLE(om_assemble assemble.cpp)

# ================
# = INSTALLATION =
# ================

install(TARGETS om_assemble om_forward om_gain om_minverser
        RUNTIME DESTINATION ${INSTALL_BIN_DIR}
        LIBRARY DESTINATION ${INSTALL_LIB_DIR}
        ARCHIVE DESTINATION ${INSTALL_LIB_DIR})

# ===========
# = TESTING =
# ===========

if (BUILD_TESTING)
    include(TestHead.cmake)

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

    # multiple nonconductive tests
    TESTHEAD(MN1)
    TESTHEAD(MN2)
    TESTHEAD(MN3)
endif()
