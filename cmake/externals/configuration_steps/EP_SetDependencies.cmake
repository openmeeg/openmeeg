################################################################################
#
# OpenMEEG
#
# Copyright (c) INRIA 2013-2017. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.
#
################################################################################

macro(EP_SetDependencies var)
    foreach(dep ${ARGN})
        if (TARGET ${dep})
            set(dependencies ${dependencies} ${dep})
        endif()
    endforeach()
    set(${var} ${dependencies})
endmacro()
