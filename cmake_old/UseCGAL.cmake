#------------------------------------------------------------
# CGAL library
#------------------------------------------------------------

option(USE_CGAL "Use CGAL meshing tools" OFF)

if (USE_CGAL)

    # what components do we want:
    set(CGAL_FIND_COMPONENTS Core)
    mark_as_advanced(CGAL_FIND_COMPONENTS)

    find_package(CGAL QUIET COMPONENTS ImageIO)
    # find_package(CGAL REQUIRED COMPONENTS Core OPTIONAL_COMPONENTS ImageIO) <- cannot since CGAL do not support OPTIONAL_COMPONENTS

    if (CGAL_ImageIO_FOUND) # for handling images (.inr format only for the moment!)
        set(CGAL_LIBRARIES CGAL::CGAL_ImageIO)
        add_definitions(-DCGAL_ImageIO)
    else()
        find_package(CGAL REQUIRED COMPONENTS ${CGAL_FIND_COMPONENTS})
        if (CGAL_FOUND)
            set(CGAL_LIBRARIES CGAL::CGAL_Core)
        else()
            message(FATAL_ERROR "Please set CGAL_DIR")
        endif()
    endif()
    set(CGAL_CXX_FLAGS ${CGAL_CXX_FLAGS_INIT} ${CGAL_SHARED_LINKER_FLAGS_INIT} ${CGAL_CXX_FLAGS_RELEASE_INIT} )
    separate_arguments(CGAL_CXX_FLAGS) # needed to remove quotes/spaces problems
    list(APPEND OpenMEEG_OTHER_LIBRARY_DIRS ${CGAL_LIBRARY_DIRS})
    list(APPEND OpenMEEG_OTHER_INCLUDE_DIRS ${CGAL_INCLUDE_DIRS})
    list(APPEND OpenMEEG_DEPENDENCIES CGAL)
    if (CGAL_3RD_PARTY_LIBRARIES)
        # old CGAL (trusty 4.2.5.ubuntu)
        set(CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_Core_LIBRARY} ${CGAL_ImageIO_LIBRARY} ${MPFR_LIBRARIES} ${GMP_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_ImageIO_3RD_PARTY_LIBRARIES})
        set(CGAL_CXX_FLAGS ${CGAL_CXX_FLAGS} ${CGAL_ImageIO_3RD_PARTY_DEFINITIONS})
    endif()
endif()
