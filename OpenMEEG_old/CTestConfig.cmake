set(CTEST_PROJECT_NAME "OpenMEEG")
set(CTEST_NIGHTLY_START_TIME "00:00:00 EST")

# Does not work currently.
# set(CTEST_DROP_METHOD "https")

set(CTEST_OUTPUT_ON_FAILURE 1)
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "cdash.inria.fr")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=openMEEG")
set(CTEST_DROP_SITE_CDASH TRUE)

set(osname ${CMAKE_SYSTEM_NAME})
set(cpu ${CMAKE_SYSTEM_PROCESSOR})
set(DISTRIB2 ${CMAKE_SYSTEM_VERSION})

set(SITE "${osname}_${DISTRIB2}_${cpu}")
set(CTEST_SITE "${osname}_${DISTRIB2}_${cpu}")
