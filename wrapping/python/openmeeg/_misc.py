# Wrappers for some function which need some parameter adjustements.

from ._openmeeg_wrapper import Head2ECoGMat_internal


def Head2ECoGMat(geom, sensors, interface_name):
    interface = geom.interface(interface_name)
    return Head2ECoGMat_internal(geom, sensors, interface)
