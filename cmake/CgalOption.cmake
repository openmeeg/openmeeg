#------------------------------------------------------------
# CGAL library
#------------------------------------------------------------

OPTION(USE_CGAL "Build the project using CGAL" OFF)
MARK_AS_ADVANCED(USE_CGAL)

IF(USE_CGAL)
    FIND_PACKAGE(CGAL QUIET COMPONENTS Core ImageIO)
    IF (CGAL_FOUND)
        INCLUDE( ${CGAL_USE_FILE} )
        SET(OPENMEEG_OTHER_LIBRARY_DIRECTORIES
            ${CGAL_LIBRARY_DIRS}
            ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES}
            )
        SET(OPENMEEG_OTHER_INCLUDE_DIRECTORIES
            ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES}
            ${CGAL_INCLUDE_DIRS}
            )
        SET(OPENMEEG_LIBRARIES
            ${OPENMEEG_LIBRARIES}
            ${CGAL_LIBRARIES}
            )
    ELSE()
        MESSAGE(FATAL_ERROR "Please set CGAL_DIR")
    ENDIF()
ENDIF()
