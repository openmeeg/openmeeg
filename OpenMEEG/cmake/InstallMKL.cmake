message(STATUS "Could not find MKL: attempting to install it from network")

#   Download MKL from Intel web site.

message(STATUS "Downloading Intel MKL 8.0, this may take a while...")

if (WIN32)
    set(MKL_URL_DIR 12079)
    set(MKL_BASE_NAME "w_mkl_2018.0.124")
    set(MKL_INSTALLER_ARCHIVE "${MKL_BASE_NAME}.exe")
elseif (APPLE)
    set(MKL_URL_DIR 12025)
    set(MKL_BASE_NAME "m_mkl_2018.0.104")
    set(MKL_INSTALLER_ARCHIVE "${MKL_BASE_NAME}.dmg")
else()
    set(MKL_URL_DIR 12070)
    set(MKL_BASE_NAME "l_mkl_2018.0.128")
    set(MKL_INSTALLER_ARCHIVE "${MKL_BASE_NAME}.tgz")
endif()

set(MKL_INSTALLER_DIR "l_mkl_2018.0.128")
set(MKL_BASE_URL "http://registrationcenter-download.intel.com/akdlm/irc_nas/tec")
file(DOWNLOAD "${MKL_BASE_URL}/${MKL_URL_DIR}/${MKL_INSTALLER_ARCHIVE}" ${CMAKE_BINARY_DIR}/${MKL_INSTALLER_ARCHIVE}
     SHOW_PROGRESS
     STATUS result)
list(GET result 0 error_code)
if (NOT ${error_code} STREQUAL "0")
    message(FATAL_ERROR "Could not download MKL install script. If no network connexion please provide MKL_DIR or environment {MKLDIR}")
endif()

#   Install MKL in a local directory.

message(STATUS "Unpacking Intel MKL")

if (WIN32)
    set(MKL_UNPACK_COMMAND ${CMAKE_BINARY_DIR}/${MKL_INSTALLER_ARCHIVE})
elseif (APPLE)
    set(MKL_UNPACK_COMMAND sudo hdiutil attach ${CMAKE_BINARY_DIR}/${MKL_INSTALLER_ARCHIVE}
                           COMMAND sudo installer -package /Volumes/${MKL_BASE_NAME}/${MKL_BASE_NAME}.pkg -target .
                           COMMAND sudo hdiutil detach /Volumes/${MKL_BASE_NAME})
else()
    # Linux
    set(MKL_UNPACK_COMMAND ${CMAKE_COMMAND} -E tar zxvf ${CMAKE_BINARY_DIR}/${MKL_INSTALLER_ARCHIVE})
endif()

message("[[execute_process(COMMAND ${MKL_UNPACK_COMMAND} OUTPUT install-mkl.out ERROR intall_mkl.err RESULT_VARIABLE mkl_unpack_result)]]")

execute_process(COMMAND ${MKL_UNPACK_COMMAND} OUTPUT_FILE ${CMAKE_BINARY_DIR}/install_mkl.out ERROR_FILE ${CMAKE_BINARY_DIR}/install_mkl.err RESULT_VARIABLE mkl_unpack_result)
file(READ ${CMAKE_BINARY_DIR}/install_mkl.out mkl_out)
message("[[mkl output: ${mkl_output}]]")
file(READ ${CMAKE_BINARY_DIR}/install_mkl.err mkl_err)
message("[[mkl err: ${mkl_err}]]")

if (NOT ${mkl_unpack_result} STREQUAL "0")
    message(FATAL_ERROR "Could not extract MKL: please look at files install-mkl.{out,err} or provide MKL_DIR or environment {MKLDIR}")
endif()

message(STATUS "Installing Intel MKL, this may take a while...")

if (WIN32)
    set(MKL_INSTALL_COMMAND "${MKL_INSTALLER_DIR}/Setup.exe install -eula=accept -output=${CMAKE_BINARY_DIR}/install-mkl.log -installdir=mkl")
else()
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/mkl)
    set(CFGFILE ${CMAKE_BINARY_DIR}/silent.cfg)
    file(APPEND ${CFGFILE} "# Generated silent configuration file\n")
    file(APPEND ${CFGFILE} "ACCEPT_EULA=accept\n")
    file(APPEND ${CFGFILE} "CONTINUE_WITH_OPTIONAL_ERROR=yes\n")
    file(APPEND ${CFGFILE} "PSET_INSTALL_DIR=${CMAKE_BINARY_DIR}/mkl\n")
    file(APPEND ${CFGFILE} "CONTINUE_WITH_INSTALLDIR_OVERWRITE=yes\n")
    file(APPEND ${CFGFILE} "COMPONENTS=ALL\n")
    file(APPEND ${CFGFILE} "PSET_MODE=install\n")
    file(APPEND ${CFGFILE} "SIGNING_ENABLED=yes\n")
    file(APPEND ${CFGFILE} "ARCH_SELECTED=ALL\n")
    set(MKL_INSTALL_COMMAND ${MKL_INSTALLER_DIR}/install.sh -s ${CFGFILE} --cli-mode --user-mode)
endif()

message("[[instal command: ${MKL_INSTALL_COMMAND}]]")

execute_process(COMMAND ${MKL_INSTALL_COMMAND}
                OUTPUT_FILE ${CMAKE_BINARY_DIR}/install-mkl.out
                ERROR_FILE ${CMAKE_BINARY_DIR}/install-mkl.err
                RESULT_VARIABLE mkl_install_result)

message("[[mkl_install_status: ${mkl_install_result}]]")
file(READ ${CMAKE_BINARY_DIR}/install-mkl.out mkl_output)
message("[[mkl output: ${mkl_output}]]")
file(READ ${CMAKE_BINARY_DIR}/install-mkl.err mkl_err)
message("[[mkl err: ${mkl_err}]]")
if (NOT ${mkl_install_result} STREQUAL "0")
    message(FATAL_ERROR "Could not install MKL: please look at files install-mkl.{out,err} or provide MKL_DIR or environment {MKLDIR}")
endif()

find_path(MKL_ROOT_DIR NAMES include/mkl_cblas.h PATHS ${CMAKE_BINARY_DIR}/mkl/mkl)
message("[[MKL_ROOT_PATH: ${MKL_ROOT_DIR}]]")
if (NOT MKL_ROOT_DIR)
    message(FATAL_ERROR "MKL seems to be incorrectly installed in ${CMAKE_BINARY_DIR}/mkl/mkl")
endif()
