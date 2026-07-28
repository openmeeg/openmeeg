macro(OpenMEEG_Executable TARGET SOURCE)
    add_executable(${TARGET} ${SOURCE})
    target_link_libraries(${TARGET} OpenMEEG::OpenMEEG)
    target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endmacro()

add_compile_options(${WERROR_COMPILE_OPTION})
OpenMEEG_Executable(om_minverser minverser.cpp)
OpenMEEG_Executable(om_gain gain.cpp)
OpenMEEG_Executable(om_forward forward.cpp)
OpenMEEG_Executable(om_assemble assemble.cpp)

# Installation

install(TARGETS om_assemble om_forward om_gain om_minverser
        ARCHIVE DESTINATION ${CMAKE_INSTALL_BINDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory(tools)
