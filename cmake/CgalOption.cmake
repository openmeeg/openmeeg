#------------------------------------------------------------
# CGAL library
#------------------------------------------------------------

option(USE_CGAL "Build the project using CGAL" OFF)
mark_as_advanced(USE_CGAL)

if (USE_CGAL)
    find_package(CGAL QUIET COMPONENTS Core ImageIO)
    if (CGAL_FOUND)
        include(${CGAL_USE_FILE} )
        set(OPENMEEG_OTHER_LIBRARY_DIRECTORIES ${CGAL_LIBRARY_DIRS} ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES})
        set(OPENMEEG_OTHER_INCLUDE_DIRECTORIES ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES} ${CGAL_INCLUDE_DIRS})
        set(OPENMEEG_LIBRARIES ${OPENMEEG_LIBRARIES} ${CGAL_LIBRARIES})
    else()
        message(FATAL_ERROR "Please set CGAL_DIR")
    endif()
endif()
