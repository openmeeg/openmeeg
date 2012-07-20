# ==============================
# = Fix DLL search pb on WIN32 =
# ==============================

IF (WIN32)
    SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${OpenMEEG_BINARY_DIR}/src")

    #   These are obsolete as per cmake 2.6

    SET(LIBRARY_OUTPUT_PATH "${OpenMEEG_BINARY_DIR}/src")
    SET(EXECUTABLE_OUTPUT_PATH "${OpenMEEG_BINARY_DIR}/src")
ENDIF()
