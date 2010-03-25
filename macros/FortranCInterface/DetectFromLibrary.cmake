#=============================================================================
# Copyright 2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

include(CheckLibraryExists)
include(CreateOutput)

function(detect_fortran_function libname library func)

#  configure_file(${FortranCInterface_SOURCE_DIR}/Input.cmake.in
#                 ${FortranCInterface_BINARY_DIR}/Input.cmake @ONLY)

#  # Detect the Fortran/C interface on the first run or when the
#  # configuration changes.
#  if(${FortranCInterface_BINARY_DIR}/Input.cmake
#      IS_NEWER_THAN ${FortranCInterface_BINARY_DIR}/Output.cmake
#      OR ${FortranCInterface_SOURCE_DIR}/Output.cmake.in
#      IS_NEWER_THAN ${FortranCInterface_BINARY_DIR}/Output.cmake
#      OR ${FortranCInterface_SOURCE_DIR}/CMakeLists.txt
#      IS_NEWER_THAN ${FortranCInterface_BINARY_DIR}/Output.cmake
#      OR ${CMAKE_CURRENT_LIST_FILE}
#      IS_NEWER_THAN ${FortranCInterface_BINARY_DIR}/Output.cmake
#      )
#    message(STATUS "Detecting Fortran/C Interface")
#  else()
#    return()
#  endif()

  # Invalidate verification results.
  unset(FortranCInterface_VERIFIED_C CACHE)
  unset(FortranCInterface_VERIFIED_CXX CACHE)

  # Detect the underscore/non underscore functions.

  string(REGEX MATCH "^.*_.*$" output ${func})
  if ("${output}" STREQUAL ${func})
    set(form "_")
  else()
    set(form "")
  endif()

  # Detect global/module functions.

  string(REGEX MATCH "^.*:.*$" output ${func})
  if ("${output}" STREQUAL ${func})
    # How to deal with with MIDDLE ??
    MESSAGE(WARNING "Modules not fully implemented yet!!")
    set(type "MODULE")
  else()
    set(type "GLOBAL")
  endif()

  string(TOUPPER "${func}" subup)
  string(TOLOWER "${func}" sublo)
  foreach(suffix "_" "")
    # foreach(base ${sublo} ${subup})
    foreach(base ${subup} ${sublo})
      foreach(prefix "" "_" "__")
        set(symbol "${prefix}${base}${suffix}")
        set(doc "function")
        message(STATUS "checking Fortran ${doc} linkage: ${symbol}")
        unset(worked CACHE)
        check_library_exists(${library} ${symbol} "" worked)
        if(worked)
          message(STATUS "found Fortran function linkage")
          set(FortranCInterface_${type}_${form}SYMBOL "${symbol}" PARENT_SCOPE)
          set(FortranCInterface_${type}_${form}PREFIX "${prefix}" PARENT_SCOPE)
          set(FortranCInterface_${type}_${form}SUFFIX "${suffix}" PARENT_SCOPE)
          if (${base} STREQUAL ${subup})
              set(FortranCInterface_${type}_${form}CASE "UPPER" PARENT_SCOPE)
          else()
              set(FortranCInterface_${type}_${form}CASE "LOWER" PARENT_SCOPE)
          endif()
          return()
        endif()
        unset(worked CACHE)
      endforeach()
    endforeach()
  endforeach()
  MESSAGE(ERROR "Function ${func} cannot be found in library ${libname}")
endfunction()
