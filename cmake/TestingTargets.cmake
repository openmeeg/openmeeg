if (BUILD_TESTING)
    add_custom_target(test
        DEPENDS build-OpenMEEG
        COMMAND cd OpenMEEG/build && make test)

    add_custom_target(test-matio
        DEPENDS build-matio
        COMMAND cd matio/build && make test)

    add_custom_target(test-all
        DEPENDS test-matio test)
endif()
