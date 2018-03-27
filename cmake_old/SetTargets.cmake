# Add Targets

set(global_targets configure install)
  
# This adds targets that will be run in each external-projects
set_property(DIRECTORY PROPERTY EP_STEP_TARGETS ${global_targets})

foreach (target ${global_targets})
    add_custom_target(${target})
endforeach()
