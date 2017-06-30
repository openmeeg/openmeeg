if (NOT HAVE_ISNORMAL_IN_NAMESPACE_STD)
    include(CheckSymbolExists)
    CHECK_SYMBOL_EXISTS(isnormal math.h HAVE_ISNORMAL_IN_MATH_H)
endif()

set(OpenMEEG_HEADERS
    analytics.h assemble.h danielsson.h domain.h forward.h gain.h geometry.h gmres.h integrator.h
    interface.h mesh.h om_utils.h om_common.h operators.h options.h PropertiesSpecialized.h vect3.h
    geometry_reader.h geometry_io.h sensors.h triangle.h Triangle_triangle_intersection.h vertex.h
#   These files are imported from another repository.
#   Please do not update them in this repository.
    DataTag.H FileExceptions.H GeometryExceptions.H Properties.H)

# install headers

install(FILES ${OpenMEEG_HEADERS}
        DESTINATION ${OpenMEEG_HEADER_INSTALLDIR} COMPONENT Development)
