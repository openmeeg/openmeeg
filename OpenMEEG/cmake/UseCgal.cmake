#------------------------------------------------------------
# CGAL library
#------------------------------------------------------------

option(USE_CGAL "Use CGAL meshing tools" OFF)

if (USE_CGAL)
    find_package(CGAL REQUIRED COMPONENTS Core ImageIO)

    if (CGAL_FOUND)
        # do not bring the whole Cgal family
        #include(${CGAL_USE_FILE})
        set(CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_Core_LIBRARY} ${CGAL_ImageIO_LIBRARY} ${MPFR_LIBRARIES} ${GMP_LIBRARIES})
        message("${CGAL_INCLUDE_DIRS} lease set ${CGAL_LIBRARY} CGAL_DIR ${CGAL_LIBRARIES}")
        set(OPENMEEG_OTHER_LIBRARY_DIRECTORIES ${CGAL_LIBRARY_DIRS} ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES})
        set(OPENMEEG_OTHER_INCLUDE_DIRECTORIES ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES} ${CGAL_INCLUDE_DIRS})
    else()
        message(FATAL_ERROR "Please set CGAL_DIR")
    endif()
endif()
