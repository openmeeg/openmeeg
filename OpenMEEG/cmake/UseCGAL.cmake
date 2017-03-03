#------------------------------------------------------------
# CGAL library
#------------------------------------------------------------

if (USE_CGAL)
    find_package(CGAL REQUIRED COMPONENTS Core ImageIO)

    if (CGAL_FOUND)
        set(CGAL_LIBRARIES CGAL::CGAL CGAL::CGAL_Core CGAL::CGAL_ImageIO)
        set(CGAL_CXX_FLAGS ${CGAL_CXX_FLAGS_INIT} ${CGAL_SHARED_LINKER_FLAGS_INIT} ${CGAL_CXX_FLAGS_RELEASE_INIT} )
        if (CGAL_3RD_PARTY_LIBRARIES)
            # old CGAL (trusty 4.2.5.ubuntu)
            set(CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_Core_LIBRARY} ${CGAL_ImageIO_LIBRARY} ${MPFR_LIBRARIES} ${GMP_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_ImageIO_3RD_PARTY_LIBRARIES})
            set(CGAL_CXX_FLAGS ${CGAL_CXX_FLAGS} ${CGAL_ImageIO_3RD_PARTY_DEFINITIONS})
        endif()
        separate_arguments(CGAL_CXX_FLAGS) # needed to remove quotes/spaces problems
        set(OPENMEEG_OTHER_LIBRARY_DIRECTORIES ${CGAL_LIBRARY_DIRS} ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES})
        set(OPENMEEG_OTHER_INCLUDE_DIRECTORIES ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES} ${CGAL_INCLUDE_DIRS})
    else()
        message(FATAL_ERROR "Please set CGAL_DIR")
    endif()
endif()
