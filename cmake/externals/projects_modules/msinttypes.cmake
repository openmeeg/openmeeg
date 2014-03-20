# OpenMEEG
#
# Copyright (c) INRIA 2013-2014. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.

function(msinttypes_project)

    set(ep msinttypes)

    # List the dependencies of the project

    set(${ep}_dependencies "")
      
    # Prepare the project

    EP_Initialisation(${ep} USE_SYSTEM OFF BUILD_SHARED_LIBS ON)

    if (NOT USE_SYSTEM_${ep})

        # Set directories

        EP_SetDirectories(${ep} EP_DIRECTORIES ep_dirs)

        # Define repository where get the sources

        if (NOT DEFINED ${ep}_SOURCE_DIR)
            set(location
                SVN_REPOSITORY "http://msinttypes.googlecode.com/svn/trunk/")
        endif()

        # Add specific cmake arguments for configuration step of the project

        set(ep_optional_args)

        # set compilation flags

        if (UNIX)
            set(${ep}_c_flags "${${ep}_c_flags} -w")
        endif()

        # Check if patch has to be applied

        ep_GeneratePatchCommand(${ep} PATCH_COMMAND)

        # Add external-project

        ExternalProject_Add(${ep}
            ${ep_dirs}
            ${location}
            ${PATCH_COMMAND}
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
        )

        # Set variable to provide infos about the project

        ExternalProject_Get_Property(${ep} binary_dir)
        set(${ep}_DIR ${binary_dir} PARENT_SCOPE)

        # Add custom targets

        EP_AddCustomTargets(${ep})

    endif()

endfunction()
