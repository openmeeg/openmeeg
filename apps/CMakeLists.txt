# macro(OPENMEEG_EXECUTABLE EXECUTABLE)
#     NEW_EXECUTABLE(${EXECUTABLE} ${ARGN} LIBRARIES OpenMEEG ${VTK_LIBRARIES} ${GIFTI_LIBRARIES})
# endmacro()

# OPENMEEG_EXECUTABLE(om_minverser minverser.cpp)
# OPENMEEG_EXECUTABLE(om_gain gain.cpp)
# OPENMEEG_EXECUTABLE(om_forward forward.cpp)
# OPENMEEG_EXECUTABLE(om_assemble assemble.cpp)

add_executable(om_minverser minverser.cpp)
target_link_libraries(om_minverser OpenMEEG::OpenMEEGMaths OpenMEEG::OpenMEEG)
target_include_directories(om_minverser PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS om_minverser
        ARCHIVE DESTINATION ${CMAKE_INSTALL_BINDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(om_gain gain.cpp)
target_link_libraries(om_gain OpenMEEG::OpenMEEGMaths OpenMEEG::OpenMEEG)
target_include_directories(om_gain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS om_gain
        ARCHIVE DESTINATION ${CMAKE_INSTALL_BINDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(om_forward forward.cpp)
target_link_libraries(om_forward OpenMEEG::OpenMEEGMaths OpenMEEG::OpenMEEG)
target_include_directories(om_forward PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS om_forward
        ARCHIVE DESTINATION ${CMAKE_INSTALL_BINDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(om_assemble assemble.cpp)
target_link_libraries(om_assemble OpenMEEG::OpenMEEGMaths OpenMEEG::OpenMEEG)
target_include_directories(om_assemble PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS om_assemble
        ARCHIVE DESTINATION ${CMAKE_INSTALL_BINDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# ================
# = INSTALLATION =
# ================

# install(TARGETS om_assemble om_forward om_gain om_minverser
#         RUNTIME DESTINATION ${INSTALL_BIN_DIR}
#         LIBRARY DESTINATION ${INSTALL_LIB_DIR}
#         ARCHIVE DESTINATION ${INSTALL_LIB_DIR})

add_subdirectory(tools)
