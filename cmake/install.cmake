# Install stuff

if (APPLE)
    set(CMAKE_MACOSX_RPATH 1)
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib/")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
elseif (UNIX)  # means if LINUX
    # to fix the installed rpath so it looks in ../lib
    # https://www.semipol.de/2012/02/16/relative-rpath-settings-with-cmake.html
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
endif()

# Find absolute path to each external lib to avoid symlinks then package it
macro(install_system_libs target)
    get_target_property(THIS_LIBS ${target} IMPORTED_LOCATION)
    foreach(LIB ${THIS_LIBS})
        get_filename_component(ABS_LIB ${LIB} REALPATH)
        list(APPEND CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${LIB} ${ABS_LIB})
    endforeach(LIB)
endmacro()

if (${MKL_USE_parallel})
    foreach(LIB ${MKL_LIBRARIES})
        if (${LIB} MATCHES "iomp")
            get_filename_component(ABS_LIB ${LIB} REALPATH)
            list(APPEND CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${LIB} ${ABS_LIB})
        endif()
    endforeach(LIB)
endif()

# install_system_libs(HDF5::HDF5)
# install_system_libs(MATIO::MATIO)

if(NOT "${EXTRA_INSTALL_LIBRARIES}" STREQUAL "")
    # https://stackoverflow.com/questions/12995166/how-can-i-install-gcc-runtime-libraries-with-cmake
    # list(APPEND CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${EXTRA_INSTALL_LIBRARIES})
    set(THIS_DESTINATION lib)
    if(WIN32)
        set(THIS_DESTINATION bin)
    endif()
    foreach(LIB ${EXTRA_INSTALL_LIBRARIES})
        message(STATUS "Adding custom ./${THIS_DESTINATION} install of ${LIB}")
        install(PROGRAMS ${LIB} DESTINATION ${THIS_DESTINATION})
    endforeach()
endif()
