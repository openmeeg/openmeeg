# quite ugly code, but it is the price for having a standalone package

# if standalone imported some libs
if (STANDALONE)
    set(imported_libs)
    list(APPEND OpenMEEG_IMPORTED_LIBS ${matio_LIBRARIES})
    list(REMOVE_DUPLICATES OpenMEEG_IMPORTED_LIBS)
    foreach(alib ${OpenMEEG_IMPORTED_LIBS})
        get_filename_component(alibname ${alib} NAME)
        string(REPLACE ${CMAKE_SHARED_LIBRARY_SUFFIX} "" alibname ${alibname})
        string(REPLACE ${CMAKE_STATIC_LIBRARY_SUFFIX} "" alibname ${alibname})
        string(REPLACE ${CMAKE_SHARED_LIBRARY_PREFIX} "" alibname ${alibname})
        # no need to export the libz nor libm nor libdl
        if (NOT ((${alibname} STREQUAL "m") OR (${alibname} STREQUAL "System.B")  OR (${alibname} STREQUAL "z")  OR (${alibname} STREQUAL "dl")))
            #add_library(${alibname} SHARED IMPORTED)
            #set_property(TARGET ${alibname} PROPERTY IMPORTED_LOCATION ${alib})
            get_filename_component(areallib ${alib} REALPATH)
            list(APPEND imported_libs ${areallib})
            if (IS_SYMLINK ${alib})
                get_filename_component(alibnameL ${alib} NAME)
                get_filename_component(areallibnameL ${areallib} NAME)
                add_custom_target(links${alibname} ALL COMMAND ln -sf ${areallibnameL} ${alibnameL})
                set(alib ${CMAKE_CURRENT_BINARY_DIR}/${alibnameL})
            endif()
            list(APPEND imported_libs ${alib})
        endif()
    endforeach()

    install(FILES ${imported_libs} DESTINATION lib
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        GROUP_EXECUTE GROUP_READ)
endif()
