# Find the matio headers and library.
#
#  matio_INCLUDE_DIRS - where to find matio.h, etc.
#  matio_LIBRARIES    - List of libraries.
#  matio_FOUND        - True if matio found.

# This module will read the variable
# MATIO_USE_STATIC_LIBRARIES to determine whether or not to prefer a
# static link to a dynamic link for MATIO and all of it's dependencies.
# To use this feature, make sure that the MATIO_USE_STATIC_LIBRARIES
# variable is set before the call to find_package.

#   We provide a module in case matio has not been found in config mode.

if (NOT matio_LIBRARIES)

    if(MATIO_USE_STATIC_LIBRARIES)
        set(HDF5_USE_STATIC_LIBRARIES TRUE)
    endif()
    find_package(HDF5 REQUIRED)

    if(MATIO_USE_STATIC_LIBRARIES AND APPLE)
        set(HDF5_LIBRARIES_XXX)
        foreach(LIB ${HDF5_LIBRARIES})
            if(${LIB} MATCHES "libsz")
                # get_filename_component(ABS_LIB ${LIB} REALPATH)
                find_library(LIBSZ
                    NAMES libsz.a
                    HINTS /usr/local/opt/szip/lib/
                    )
                set(LIB ${LIBSZ})
            endif()
            set(HDF5_LIBRARIES_XXX ${HDF5_LIBRARIES_XXX} ${LIB})
        endforeach(LIB)
        set(HDF5_LIBRARIES ${HDF5_LIBRARIES_XXX})
    endif()

    # Make a modern cmake interface to HDF5
    add_library(HDF5::HDF5 INTERFACE IMPORTED)
    set_target_properties(HDF5::HDF5 PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES}")

    # Look for the header file.
    set(conda_matio /home/travis/miniconda/pkgs/libmatio-1.5.12-0/)
    find_path(matio_INCLUDE_DIR
	    HINTS
        	$ENV{matio_dir}include
          ${conda_matio}include
	    NAMES
	    	matio.h
	    )

    message(STATUS "matio.h ${matio_INCLUDE_DIR}")
    mark_as_advanced(matio_INCLUDE_DIR)

    # Look for the library.

    # XXXX This needs to go out !
    set(matio_LIB_SEARCH_PATHS
        C:/conda/Library/
        C:/conda/Library/lib
        C:/conda/Library/bin
        $ENV{matio_dir}
        $ENV{matio_dir}lib
        $ENV{matio_dir}bin
        ${conda_matio}
        ${conda_matio}lib
        ${conda_matio}bin
        )

    set(MATIO_NAMES matio libmatio)
    if(MATIO_USE_STATIC_LIBRARIES)
        set(MATIO_NAMES  libmatio.a ${MATIO_NAMES})
    endif()

    find_library(matio_LIBRARY
        HINTS
            ${matio_LIB_SEARCH_PATHS}
        NAMES
            ${MATIO_NAMES}
        )
    message(STATUS "matio_library ${matio_LIBRARY}")
    mark_as_advanced(matio_LIBRARY)

    # handle the QUIETLY and REQUIRED arguments and set matio_FOUND to TRUE if
    # all listed variables are TRUE
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(matio DEFAULT_MSG matio_LIBRARY matio_INCLUDE_DIR)

    if (matio_FOUND AND NOT TARGET MATIO::MATIO)
        add_library(MATIO::MATIO UNKNOWN IMPORTED)
        set_target_properties(MATIO::MATIO PROPERTIES
            IMPORTED_LINK_LIBRARIES HDF5::HDF5
            INTERFACE_INCLUDE_DIRECTORIES ${matio_INCLUDE_DIR}
            IMPORTED_LOCATION ${matio_LIBRARY}
        )
    endif()
endif()
