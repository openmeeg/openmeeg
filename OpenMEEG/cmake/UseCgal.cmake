#------------------------------------------------------------
# CGAL library
#------------------------------------------------------------

option(USE_CGAL "Use CGAL meshing tools" OFF)

if (USE_CGAL)
    find_package(CGAL REQUIRED COMPONENTS Core ImageIO)

    if (CGAL_FOUND)
        # do not bring the whole Cgal family
        set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE)
        set(CGAL_HEADER_ONLY TRUE)
        include(${CGAL_USE_FILE})
        set(OPENMEEG_OTHER_LIBRARY_DIRECTORIES ${CGAL_LIBRARY_DIRS} ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES})
        set(OPENMEEG_OTHER_INCLUDE_DIRECTORIES ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES} ${CGAL_INCLUDE_DIRS})
        # no need for this
        # set(OPENMEEG_LIBRARIES ${OPENMEEG_LIBRARIES} ${CGAL_LIBRARIES})
    else()
        message(FATAL_ERROR "Please set CGAL_DIR")
    endif()
endif()
