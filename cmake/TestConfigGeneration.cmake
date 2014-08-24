include(GetPrerequisites)

get_prerequisites("${OM_ASSEMBLE_EXE}" PREREQUISITES 1 0 "${OM_ASSEMBLE_EXE}" "")
foreach (i ${PREREQUISITES})
    get_filename_component(realpath ${i} REALPATH)
    get_filename_component(dir ${realpath} PATH)
    set(LIBRARY_PATHS "${dir};${LIBRARY_PATHS}")
endforeach()
file(WRITE ${CMAKE_BINARY_DIR}/TestConfig.cmake "set(ENV{PATH} \"${LIBRARY_PATHS}\$ENV{PATH}\")\n")
