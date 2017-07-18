#------------------------------------------------------------
# OpenBLAS library
#------------------------------------------------------------

if (USE_OPENBLAS)
    if (WIN32)
        # special case for windows
        if (${CMAKE_PROJECT_NAME} STREQUAL "OpenMEEG")
            # so that the dlls can be found at runtime
            set(CMAKE_MSVCIDE_RUN_PATH ${OpenBLAS_DIR}/bin ${OpenBLAS_DIR}/../mingw32_dll ${OpenBLAS_DIR}/../mingw64_dll ${CMAKE_MSVCIDE_RUN_PATH})
            # add openblas as a shared imported lib
            add_library(libopenblas SHARED IMPORTED)
            set_property(TARGET libopenblas PROPERTY IMPORTED_LOCATION ${OpenBLAS_DIR}/bin/libopenblas.dll)
            set_property(TARGET libopenblas PROPERTY IMPORTED_IMPLIB   ${OpenBLAS_DIR}/lib/libopenblas.dll.a)

            file(GLOB MinGW_LIBS "${OpenBLAS_DIR}/../mingw*/lib*")

            install(FILES ${MinGW_LIBS} ${OpenBLAS_DIR}/bin/libopenblas.dll DESTINATION bin
                    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_EXECUTE GROUP_READ)
        endif()
        # set variables manually as the OpenBLASConfig.cmake given is wrong
        set(OpenBLAS_INCLUDE_DIR ${OpenBLAS_DIR}/include)
        set(OpenBLAS_LIBRARIES libopenblas)
        set(OpenBLAS_FOUND TRUE)
    else()
        find_package(OpenBLAS ${FIND_MODE} MODULE)
        # if we are in a standalone mode
        if ((${CMAKE_PROJECT_NAME} STREQUAL "OpenMEEG") AND STANDALONE)
            # add openblas and lapacke as a shared imported lib
            add_library(libopenblas SHARED IMPORTED)
            set_property(TARGET libopenblas PROPERTY IMPORTED_LOCATION ${OpenBLAS_LIB})
            if (Lapacke_LIB)
                add_library(lapacke SHARED IMPORTED)
                set_property(TARGET lapacke PROPERTY IMPORTED_LOCATION ${Lapacke_LIB})
            endif()

            install(FILES ${OpenBLAS_LIBRARIES} DESTINATION lib
                    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_EXECUTE GROUP_READ)
        endif()
    endif()
    if (OpenBLAS_FOUND)
        include_directories(${OpenBLAS_INCLUDE_DIR})
        set(LAPACK_LIBRARIES ${OpenBLAS_LIBRARIES})
        list(APPEND OpenMEEG_DEPENDENCIES OpenBLAS)
    endif()
endif()
