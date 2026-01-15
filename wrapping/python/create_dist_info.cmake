set(OPENMEEG_VERSION $ENV{SETUPTOOLS_SCM_PRETEND_VERSION})
message("${OPENMEEG_VERSION}")
message("$ENV{BDIR}")
message("Python: Generating $ENV{BDIR}/openmeeg-${OPENMEEG_VERSION}.dist_info")
execute_process(WORKING_DIRECTORY $ENV{BDIR}
                COMMAND $ENV{PYTHON} $ENV{SDIR}/setup.py dist_info --output-dir $ENV{BDIR})
