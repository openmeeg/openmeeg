# quite ugly code, but it is the price for having a standalone package

# if standalone imported some libs
if (STANDALONE)
    set(imported_libs)
    list(APPEND OpenMEEG_IMPORTED_LIBS ${matio_LIBRARIES})
    list(REMOVE_DUPLICATES OpenMEEG_IMPORTED_LIBS)
    foreach(alib ${OpenMEEG_IMPORTED_LIBS})
        get_filename_component(alibname ${alib} NAME)
        get_filename_component(areallib ${alib} REALPATH)
        get_filename_component(areallibname ${areallib} NAME)
        string(REPLACE ${CMAKE_SHARED_LIBRARY_SUFFIX} "" alibname ${alibname})
        string(REPLACE ${CMAKE_STATIC_LIBRARY_SUFFIX} "" alibname ${alibname})
        string(REPLACE ${CMAKE_SHARED_LIBRARY_PREFIX} "" alibname ${alibname})
        # no need to export the libz nor libm nor libdl nor System.B
        if (NOT ((${alibname} STREQUAL "m") OR (${areallibname} MATCHES "System.B")
            OR (${alibname} STREQUAL "z")  OR (${alibname} STREQUAL "dl")))
            #add_library(${alibname} SHARED IMPORTED) # not needed as too many libs
            #set_property(TARGET ${alibname} PROPERTY IMPORTED_LOCATION ${alib})
            list(APPEND imported_libs ${areallib})
            if (IS_SYMLINK ${alib})
                message("alib: ${alib}")
                message("alibname: ${alibname}")
                message("areallib: ${areallib}")
                message("areallibname: ${areallibname}")
                get_filename_component(alibnameL ${alib} NAME)
                # in case we need to do an extra sym link from libmatio.so.2 -> libmatio.so.2.0.2
                get_filename_component(areallibnameExt ${areallib} EXT)
                string(REGEX MATCHALL "\\.[^.]*" exts ${areallibnameExt})
                list(LENGTH exts nb_ext)
                if (${nb_ext} GREATER 2)
                    list(GET exts 1 outt)
                    set(linkname ${alibnameL}${outt})
                    if (APPLE)
                        list(GET exts 0 outt)
                        set(linkname ${CMAKE_SHARED_LIBRARY_PREFIX}${alibname}${outt}${CMAKE_SHARED_LIBRARY_SUFFIX})
                    endif()
                    add_custom_target(links${linkname} ALL
                        COMMAND ln -sf ${areallibname} ${linkname}
                        COMMENT "linking ${linkname} -> ${areallibname}")
                    list(APPEND imported_libs ${CMAKE_CURRENT_BINARY_DIR}/${linkname})
                endif()
                # sym link from libmatio.so -> libmatio.so.2.0.2
                add_custom_target(links${alibnameL} ALL
                    COMMAND ln -sf ${areallibname} ${alibnameL}
                    COMMENT "linking ${alibnameL} -> ${areallibname}")
                list(APPEND imported_libs ${CMAKE_CURRENT_BINARY_DIR}/${alibnameL})
            endif()
        endif()
    endforeach()

    install(FILES ${imported_libs} DESTINATION lib
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        GROUP_EXECUTE GROUP_READ)
endif()
