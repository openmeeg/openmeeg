macro(SetRevision major minor patch dev)

    set(${PROJECT_NAME}_VERSION_MAJOR ${major} CACHE STRING "OpenMEEG-superproject major version number.")
    mark_as_advanced(${PROJECT_NAME}_VERSION_MAJOR)

    set(${PROJECT_NAME}_VERSION_MINOR ${minor} CACHE STRING "OpenMEEG-superproject minor version number.")
    mark_as_advanced(${PROJECT_NAME}_VERSION_MINOR)

    set(${PROJECT_NAME}_VERSION_PATCH ${patch} CACHE STRING "OpenMEEG-superproject build version number.")
    mark_as_advanced(${PROJECT_NAME}_VERSION_PATCH)

    set(${PROJECT_NAME}_VERSION_TWEAK ${dev} CACHE STRING "OpenMEEG-superproject development marker.")
    mark_as_advanced(${PROJECT_NAME}_VERSION_TWEAK)

    if (NOT ${${PROJECT_NAME}_VERSION_TWEAK} STREQUAL "")
        set(${PROJECT_NAME}_VERSION 
            ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}.${${PROJECT_NAME}_VERSION_TWEAK})
    else()
        set(${PROJECT_NAME}_VERSION
            ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH})
    endif()
endmacro()
