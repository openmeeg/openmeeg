#------------------------------------------------------------
# VTK library
#------------------------------------------------------------

option(USE_VTK "Build the project using VTK" OFF)

if (USE_VTK)
    find_package(VTK)
    if (VTK_FOUND)
        include(${VTK_USE_FILE})
        set (VTK_LIBRARIES vtkGraphics vtkIO vtkCommon)
        # vtkRendering vtkHybrid vtkImaging vtkFiltering vtkGenericFiltering vtkDICOMParser vtkzlib
        set(OPENMEEG_OTHER_LIBRARY_DIRECTORIES ${VTK_LIBRARY_DIRS} ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES})
        set(OPENMEEG_OTHER_INCLUDE_DIRECTORIES ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES} ${VTK_INCLUDE_DIRS})
        set(OPENMEEG_LIBRARIES ${OPENMEEG_LIBRARIES} ${VTK_LIBRARIES})
    else()
        message(FATAL_ERROR "Please set VTK_DIR")
    endif()
endif()
