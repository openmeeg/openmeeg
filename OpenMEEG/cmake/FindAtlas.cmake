
set(ATLAS_LIB_SEARCHPATH
    /usr/lib64/atlas/sse3 /usr/lib/atlas/sse3 /usr/lib/sse3
    /usr/lib64/atlas/sse2 /usr/lib/atlas/sse2 /usr/lib/sse2
    /usr/lib64/atlas /usr/lib/atlas /usr/lib/atlas-base /usr/lib64/atlas-base
    /usr/lib64/ /usr/lib/)

macro(find_atlas_lib)
    set(lib_found FALSE)
    foreach (LIB ${ARGN})
        message("Searching: ${LIB}")
        set(LIBNAMES ${LIB})
        if (${LIB} STREQUAL "clapack")
            set(LIBNAMES ${LIB} lapack_atlas)
        endif()
        find_library(${LIB}_PATH
            NAMES ${LIBNAMES}
            PATHS ${ATLAS_LIB_SEARCHPATH}
            NO_DEFAULT_PATH
            NO_CMAKE_ENVIRONMENT_PATH
            NO_CMAKE_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH)
        if (${LIB}_PATH)
            get_filename_component(LAPACK_ROOT_DIR ${${LIB}_PATH} DIRECTORY)
            set(Atlas_LIBRARIES ${Atlas_LIBRARIES} ${${LIB}_PATH})
            mark_as_advanced(${LIB}_PATH)
            set(lib_found TRUE)
            break()
        endif()
    endforeach()
endmacro()

find_path(Atlas_INCLUDE_DIR clapack.h /usr/include/atlas /usr/include/ NO_DEFAULT_PATH)

find_atlas_lib(tatlas satlas atlas)
find_atlas_lib(clapack)
find_atlas_lib(lapack)
find_atlas_lib(cblas f77blas blas)

if(Atlas_LIBRARIES AND Atlas_INCLUDE_DIR)
    set(Atlas_FOUND true)
    mark_as_advanced(Atlas_INCLUDE_DIR)
else()
    unset(Atlas_INCLUDE_DIR CACHE)
endif()
