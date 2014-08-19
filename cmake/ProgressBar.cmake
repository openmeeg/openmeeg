#-----------------------------------------------
# Progress bar
#-----------------------------------------------

option(USE_PROGRESSBAR "Show ascii progress bar when assembling matrices" OFF)

if (USE_PROGRESSBAR)
    add_definitions(-DUSE_PROGRESSBAR)
endif()
