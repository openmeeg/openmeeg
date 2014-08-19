if (MSINTTYPES)
    add_custom_target(uninstall-msinttypes
        COMMAND cd msinttypes/build && make uninstall)
endif()

add_custom_target(uninstall-zlib
    COMMAND cd zlib/build && make uninstall)

if (LAPACK)
    add_custom_target(uninstall-lapack
        COMMAND cd clapack/build && make uninstall)
endif()

add_custom_target(uninstall-hdf5
    COMMAND cd hdf5/build && make uninstall)

add_custom_target(uninstall-matio
    COMMAND cd matio/build && make uninstall)

add_custom_target(uninstall-OpenMEEG
    COMMAND cd OpenMEEG/build && make uninstall)

add_custom_target(uninstall
    DEPENDS uninstall-hdf5 uninstall-matio uninstall-OpenMEEG)
