################################################################################
#
# medInria
#
# Copyright (c) INRIA 2013. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.
#
################################################################################

function(ep_SetDirectories ep  
  EP_DIRECTORIES ep_dirs
  )

## #############################################################################
## Define a directory for each target of the project
## #############################################################################

set(DIR_VAR_NAMES 
  DOWNLOAD 
  BINARY 
  STAMP 
  INSTALL 
  TMP
  )

set(DIR_NAMES     
  ""       
  build  
  stamp 
  install/${CMAKE_CFG_INTDIR} 
  tmp
  )

set(dirs PREFIX ${ep})
foreach(i RANGE 4)
  list(GET DIR_VAR_NAMES ${i} var)
  list(GET DIR_NAMES     ${i} dir)
  set(dirs ${dirs} ${var}_DIR ${ep}/${dir})
endforeach()

## #############################################################################
## Look for and define the source directory of the project 
## #############################################################################

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ep}/CMakeLists.txt 
OR EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ep}/configure)
  set(${ep}_SOURCE_DIR SOURCE_DIR ${CMAKE_SOURCE_DIR}/${ep} PARENT_SCOPE)    
endif()

set(source_dir ${CMAKE_SOURCE_DIR}/${ep})
set(source_dir ${source_dir} PARENT_SCOPE) 
set(dirs ${dirs} SOURCE_DIR ${source_dir})
set(${ep_dirs} ${dirs} PARENT_SCOPE) 

endfunction()
