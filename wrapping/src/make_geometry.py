# Build a geometry with given interfaces and domains.


def make_geometry(meshes, interfaces, domains):
    """Make a geometry from dictionary of domains and a list of interfaces


    Parameters
    ----------
    interfaces : list
        The XXX
    domains : dict
        The domains.

    Returns
    -------
    geometry : isinstance of om.Geometry
        The geometry that can be used in OpenMEEG.
    """

    if not isinstance(meshes, dict) or len(meshes) == 0:
        raise ValueError(
            "Wrong argument (should be a non empty dictionary of named meshes)."
            " Got %s" % type(meshes))

    if not isinstance(interfaces, dict) or len(interfaces) == 0:
        raise ValueError(
            "Wrong argument (should be a non empty dictionary of named interfaces)"
            " Got %s" % type(interfaces))

    if not isinstance(domains, dict) or len(domains) == 0:
        raise ValueError(
            "Wrong argument (should be a non empty dictionary of named domains)"
            " Got %s" % type(domains))

    # First add mesh points

    indmaps = dict()
    geom = Geometry()  # noqa
    for name, mesh in meshes.items():
        print(name, mesh)
        indmaps[name] = geom.add_vertices(mesh[0])

    # Create meshes

    for name, mesh in meshes.items():
        om_mesh = geom.add_mesh(name)
        om_mesh.add_triangles(mesh[1], indmaps[name])
        om_mesh.update(True)

    for dname, domain in domains.items():
        domain_interfaces, conductivity = domain

        if not isinstance(domain_interfaces, list) or len(domain_interfaces) == 0:
            raise Exception("wrong description of domain (" +
                            dname + "), should be a non-empty list of interfaces")

        om_domain = Domain(dname)  # noqa
        om_domain.set_conductivity(conductivity)

        for iname, side in domain_interfaces:
            if iname not in interfaces:
                raise Exception(
                    "Domain " + dname + " contains and unknown interface " + iname + ".")
            oriented_meshes = interfaces[iname]
            if type(oriented_meshes) != list or len(oriented_meshes) == 0:
                raise Exception("Interface definition " + iname +
                                " first argument should be a non-empty list of (mesh,orientation)")
            if side != SimpleDomain.Inside and side != SimpleDomain.Outside:  # noqa
                raise Exception(
                    "Domain " + dname + ": interface " + iname + " has a wrong side direction (In/Out)")

            om_interface = Interface(iname)  # noqa
            for mesh_name, orientation in oriented_meshes:
                if orientation != OrientedMesh.Normal and orientation != OrientedMesh.Opposite:  # noqa
                    raise Exception(
                        "Wrong description for interface (" + iname + "), second tuple member should a be an orientation")

                mesh = geom.mesh(mesh_name)
                oriented_mesh = OrientedMesh(mesh, orientation)  # noqa
                om_interface.oriented_meshes().push_back(oriented_mesh)  # noqa

            om_domain.boundaries().push_back(SimpleDomain(om_interface, side))  # noqa
        geom.domains().push_back(om_domain)

    for dname, domain in domains.items():
        domain_interfaces, conductivity = domain
        om_domain = Domain(dname)  # noqa
        om_domain.set_conductivity(conductivity)
        for iname, side in domain_interfaces:
            oriented_meshes = interfaces[iname]
            om_interface = Interface(iname)  # noqa
            for mesh, orientation in oriented_meshes:
                om_mesh = geom.mesh(mesh)
                oriented_mesh = OrientedMesh(om_mesh, orientation)  # noqa
                om_interface.oriented_meshes().push_back(oriented_mesh)
            om_domain.boundaries().push_back(SimpleDomain(om_interface, side))  # noqa
        geom.domains().push_back(om_domain)

    geom.finalize()
    return geom
