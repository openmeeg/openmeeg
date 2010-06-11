# ==============================
# = Fix DLL search pb on WIN32 =
# ==============================

IF (WIN32)
    SET(LIBRARY_OUTPUT_PATH "${OpenMEEG_BINARY_DIR}/src")
    SET(EXECUTABLE_OUTPUT_PATH "${OpenMEEG_BINARY_DIR}/src")
ENDIF()
