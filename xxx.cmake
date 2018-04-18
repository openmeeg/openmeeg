
####
# MATIO-Windows shipping work-around
# THIS NEEDS TO GO OUT
#####
if(WIN32)
  set(xxxx "C:/conda/Library/bin/libmatio.dll")

  # find_path(xxxx
  #     HINTS
  #         ${xxxx_LIB_SEARCH_PATHS}
  #     NAMES
  #         libmatio
  #     NO_DEFAULT_PATH
  #     )
  message(DEBUG xxxx: ${xxxx})
  if(EXISTS ${xxxx})
    install(FILES ${xxxx} DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
  endif()
endif(WIN32)
