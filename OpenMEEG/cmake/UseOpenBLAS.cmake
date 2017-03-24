#------------------------------------------------------------
# OpenBLAS library
#------------------------------------------------------------

if (USE_OPENBLAS)
    if(WIN32)
        if (NOT ${CMAKE_PROJECT_NAME} STREQUAL "OpenMEEG")
            # download Openblas
            if (FORCE_BUILD_32BITS)
                if (NOT EXISTS ${CMAKE_SOURCE_DIR}/OpenBLAS.zip)
                    file(DOWNLOAD "https://downloads.sourceforge.net/project/openblas/v0.2.19/OpenBLAS-v0.2.19-Win32.zip" ${CMAKE_SOURCE_DIR}/OpenBLAS.zip SHOW_PROGRESS EXPECTED_MD5 "cc29b41bc1fe41e8ef5ecb452e1cd70c")
                    # file(DOWNLOAD "https://downloads.sourceforge.net/project/openblas/v0.2.14/OpenBLAS-v0.2.14-Win32.zip" ${CMAKE_SOURCE_DIR}/OpenBLAS.zip SHOW_PROGRESS EXPECTED_MD5 "eefdf170439620d78fabb3139b7aeb2f")
                endif()
                if (NOT EXISTS ${CMAKE_SOURCE_DIR}/mingw.zip)
                    file(DOWNLOAD "https://sourceforge.net/projects/openblas/files/v0.2.14/mingw32_dll.zip" ${CMAKE_SOURCE_DIR}/mingw.zip SHOW_PROGRESS EXPECTED_MD5 "47f8b18b7b99ea0ca452bbfc4f6ef579")
                endif()
            else()
                if (NOT EXISTS ${CMAKE_SOURCE_DIR}/OpenBLAS.zip)
                    file(DOWNLOAD "https://downloads.sourceforge.net/project/openblas/v0.2.19/OpenBLAS-v0.2.19-Win64-int32.zip" ${CMAKE_SOURCE_DIR}/OpenBLAS.zip SHOW_PROGRESS EXPECTED_MD5 "7ff6092397a93494c137e23670dd72ec")
                    # file(DOWNLOAD "https://downloads.sourceforge.net/project/openblas/v0.2.14/OpenBLAS-v0.2.19-Win64-int32.zip" ${CMAKE_SOURCE_DIR}/OpenBLAS.zip SHOW_PROGRESS EXPECTED_MD5 "bb59507959975d8d55f3e7eb0ecd5ea3")
                endif()
                if (NOT EXISTS ${CMAKE_SOURCE_DIR}/mingw.zip)
                    file(DOWNLOAD "https://sourceforge.net/projects/openblas/files/v0.2.14/mingw64_dll.zip" ${CMAKE_SOURCE_DIR}/mingw.zip SHOW_PROGRESS EXPECTED_MD5 "e619f1f936638240472397e0f7970e66")
                endif()
            endif()
            # extract the openblas/mingw archives
            execute_process(COMMAND ${CMAKE_COMMAND} -E tar xfz ${CMAKE_SOURCE_DIR}/OpenBLAS.zip)
            execute_process(COMMAND ${CMAKE_COMMAND} -E tar xfz ${CMAKE_SOURCE_DIR}/mingw.zip)
        else()
            file(GLOB OpenBLAS_DIR "${CMAKE_BINARY_DIR}/../../OpenBLAS*")
            file(GLOB MinGW_LIBS "${CMAKE_BINARY_DIR}/../../mingw*/lib*")

            # add openblas as a shared imported lib
            add_library(libopenblas SHARED IMPORTED)
            set_property(TARGET libopenblas PROPERTY IMPORTED_LOCATION ${OpenBLAS_DIR}/bin/libopenblas.dll)
            set_property(TARGET libopenblas PROPERTY IMPORTED_IMPLIB   ${OpenBLAS_DIR}/lib/libopenblas.dll.a)

            # copy Dlls to build and install.
            set(OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/src/${CMAKE_BUILD_TYPE})
            foreach(lib ${MinGW_LIBS} ${OpenBLAS_DIR}/bin/libopenblas.dll)
               configure_file(${lib} ${OUTPUT_DIRECTORY} COPYONLY)
            endforeach()

            install(FILES ${MinGW_LIBS} ${OpenBLAS_DIR}/bin/libopenblas.dll DESTINATION bin
                    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_EXECUTE GROUP_READ)

            # set variables manually as the OpenBLASConfig.cmake given is wrong
            set(OpenBLAS_INCLUDE_DIR ${OpenBLAS_DIR}/include)
        endif()
        set(OpenBLAS_LIBRARIES libopenblas)
        set(OpenBLAS_FOUND TRUE)
    else()
        find_package(OpenBLAS ${REQUIRED} MODULE)
    endif()
    if (OpenBLAS_FOUND)
        include_directories(${OpenBLAS_INCLUDE_DIR})
        set(LAPACK_LIBRARIES ${OpenBLAS_LIBRARIES})
    endif()
endif()
