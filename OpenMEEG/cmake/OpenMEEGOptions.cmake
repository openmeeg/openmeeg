include(BlasLapackOption)

option(ENABLE_PACKAGING "Enable Packaging" OFF)
option(ENABLE_PYTHON "Enable Python Wrapping" ON)
option(USE_OMP "Use OpenMP" OFF)
option(USE_GIFTI "Use GIFTI IO support" OFF)
option(USE_VTK "Use VTK" OFF)
option(USE_CGAL "Use CGAL meshing tools" OFF)
option(USE_KWSTYLE "Checking code syntax using KWStyle" OFF)
option(BUILD_TESTING "Build the testing tree" ON)
option(BUILD_DOCUMENTATION "Build the documentation" ON)

mark_as_advanced(USE_KWSTYLE)
mark_as_advanced(USE_GIFTI)
