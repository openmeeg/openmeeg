#------------------------------------------------------------
# VTK library
#------------------------------------------------------------

OPTION(USE_VTK "Build the project using VTK" OFF)

IF (USE_VTK)
    FIND_PACKAGE(VTK)
    IF (VTK_FOUND)
        INCLUDE(${VTK_USE_FILE})
        SET (VTK_LIBRARIES
            # vtkRendering
            vtkGraphics
            # vtkHybrid
            # vtkImaging
            vtkIO
            # vtkFiltering
            # vtkGenericFiltering
            vtkCommon
            # vtkDICOMParser
            # vtkzlib
        )
        SET(OPENMEEG_OTHER_LIBRARY_DIRECTORIES
            ${VTK_LIBRARY_DIRS}
            ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES}
        )
        SET(OPENMEEG_OTHER_INCLUDE_DIRECTORIES
            ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES}
            ${VTK_INCLUDE_DIRS}
        )
        SET(OPENMEEG_LIBRARIES
            ${OPENMEEG_LIBRARIES}
            ${VTK_LIBRARIES}
        )
    ELSE()
        MESSAGE(FATAL_ERROR "Please set VTK_DIR")
    ENDIF()
ENDIF()
