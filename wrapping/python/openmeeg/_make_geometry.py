# Build a geometry with given interfaces and domains.
from pathlib import Path

import numpy as np

from ._openmeeg_wrapper import (
    Domain,
    Geometry,
    Interface,
    Mesh,
    OrientedMesh,
    SimpleDomain,
)


def _mesh_vertices_and_triangles(mesh):
    mesh_vertices = mesh.geometry().vertices()
    vertices = np.array([vertex.array() for vertex in mesh_vertices], np.float64)
    mesh_triangles = mesh.triangles()
    triangles = np.array(
        [mesh.triangle(triangle).array() for triangle in mesh_triangles],
        dtype=np.int64,
    )
    return vertices, triangles


def make_geometry(meshes, interfaces, domains):
    """Make a geometry from dictionary of domains and a list of interfaces.

    Parameters
    ----------
    meshes : dict
        Dictionary of meshes, indexed by domain name. Meshes can be
        either instances of Mesh or tuples of (vertices, triangles).
    interfaces : dict
        Dictionary of interfaces, indexed by interface name.
    domains : dict
        The domains.

    Returns
    -------
    geometry : isinstance of om.Geometry
        The geometry that can be used in OpenMEEG.
    """
    if not isinstance(meshes, dict) or len(meshes) == 0:
        raise ValueError(
            "Wrong argument (should be a non empty dictionary of named "
            f"meshes). Got {type(meshes)}"
        )

    if not isinstance(interfaces, dict) or len(interfaces) == 0:
        raise ValueError(
            "Wrong argument (should be a non empty dictionary of named "
            f"interfaces). Got {type(interfaces)}"
        )

    if not isinstance(domains, dict) or len(domains) == 0:
        raise ValueError(
            "Wrong argument (should be a non empty dictionary of named "
            f"domains). Got {type(domains)}"
        )

    # Normalize mesh inputs to numpy arrays

    for name, mesh in meshes.items():
        if isinstance(mesh, Mesh):
            meshes[name] = _mesh_vertices_and_triangles(mesh)
        elif isinstance(mesh, (list, tuple)):
            pass
        else:
            raise ValueError(
                f"Wrong argument (should be a Mesh or a tuple of "
                f"vertices and triangles). Got {type(mesh)}"
            )

    # First add mesh points

    indmaps = dict()
    geom = Geometry(len(meshes))
    for name, mesh in meshes.items():
        indmaps[name] = geom.add_vertices(mesh[0])

    # Create meshes

    for name, mesh in meshes.items():
        om_mesh = geom.add_mesh(name)
        om_mesh.add_triangles(mesh[1], indmaps[name])
        om_mesh.update(True)

    del meshes, indmaps

    for dname, domain in domains.items():
        domain_interfaces, conductivity = domain

        if not isinstance(domain_interfaces, list) or len(domain_interfaces) == 0:
            raise Exception(
                f"wrong description of domain ({dname}), should be a "
                "non-empty list of interfaces"
            )

        om_domain = Domain(dname)
        om_domain.set_conductivity(conductivity)

        for iname, side in domain_interfaces:
            if iname not in interfaces:
                raise Exception(
                    f"Domain {dname} contains and unknown interface {iname}."
                )
            oriented_meshes = interfaces[iname]
            if not isinstance(oriented_meshes, list) or len(oriented_meshes) == 0:
                raise Exception(
                    f"Interface definition {iname} first argument should be a "
                    "non-empty list of (mesh,orientation)"
                )
            if side != SimpleDomain.Inside and side != SimpleDomain.Outside:
                raise Exception(
                    f"Domain {dname}: interface {iname} has a wrong side "
                    "direction (In/Out)"
                )

            om_interface = Interface(iname)
            for mesh_name, orientation in oriented_meshes:
                if (
                    orientation != OrientedMesh.Normal
                    and orientation != OrientedMesh.Opposite
                ):
                    raise Exception(
                        f"Wrong description for interface ({iname}), second "
                        "tuple member should a be an orientation"
                    )

                mesh = geom.mesh(mesh_name)
                oriented_mesh = OrientedMesh(mesh, orientation)
                om_interface.oriented_meshes().push_back(oriented_mesh)

            om_domain.boundaries().push_back(SimpleDomain(om_interface, side))
        geom.domains().push_back(om_domain)

    geom.finalize()
    return geom


def make_nested_geometry(meshes, conductivity):
    """Make a geometry from a list of meshes assumed to be nested.

    Parameters
    ----------
    meshes : list
        List of meshes from inner to outer. For now only 1 or 3 layer
        models are supported.
    conductivity : array-like
        The list of conductivities for each domain from the inside to the
        outside. For example [1, 0.0125, 1].

    Returns
    -------
    geometry : isinstance of om.Geometry
        The geometry that can be used in OpenMEEG.
    """
    if not isinstance(meshes, list) or len(meshes) not in [1, 3]:
        raise ValueError(
            "Wrong argument (should be a list of 1 or 3 meshes). " f"Got {type(meshes)}"
        )

    if len(meshes) == 3:
        # Convert meshes to dictionary of meshes for make_geometry
        meshes = {
            "Cortex": meshes[0],
            "Skull": meshes[1],
            "Head": meshes[2],
        }
        brain_conductivity, skull_conductivity, scalp_conductivity = conductivity

        # It should be possible to have multiple oriented meshes per interface.
        # e.g.
        # interface1 = [(m1,om.OrientedMesh.Normal),
        #               (m2,om.OrientedMesh.Opposite),
        #               (m3,om.OrientedMesh.Normal)]
        # It should also be possible to have a name added at the beginning of the
        # tuple.

        interfaces = {
            "Cortex": [("Cortex", OrientedMesh.Normal)],
            "Skull": [("Skull", OrientedMesh.Normal)],
            "Head": [("Head", OrientedMesh.Normal)],
        }

        domains = {
            "Scalp": (
                [
                    ("Skull", SimpleDomain.Outside),
                    ("Head", SimpleDomain.Inside),
                ],
                scalp_conductivity,
            ),
            "Brain": ([("Cortex", SimpleDomain.Inside)], brain_conductivity),
            "Air": ([("Head", SimpleDomain.Outside)], 0.0),
            "Skull": (
                [
                    ("Cortex", SimpleDomain.Outside),
                    ("Skull", SimpleDomain.Inside),
                ],
                skull_conductivity,
            ),
        }
    else:  # assume 1 layer model
        # Convert meshes to dictionary of meshes for make_geometry
        meshes = {
            "Cortex": meshes[0],
        }
        (brain_conductivity,) = conductivity

        interfaces = {
            "Cortex": [("Cortex", OrientedMesh.Normal)],
        }

        domains = {
            "Brain": ([("Cortex", SimpleDomain.Inside)], brain_conductivity),
            "Air": ([("Cortex", SimpleDomain.Outside)], 0.0),
        }

    geom = make_geometry(meshes, interfaces, domains)
    return geom


def read_geometry(geom_file, cond_file):
    """Read a geometry from a file.

    Parameters
    ----------
    geom_file : path-like
        The name of the geometry file.
    cond_file : path-like
        The name of the conductivity file.

    Returns
    -------
    geometry : isinstance of om.Geometry
        The geometry that can be used in OpenMEEG.
    """
    geom_file = str(Path(geom_file).resolve())
    cond_file = str(Path(cond_file).resolve())
    return Geometry(geom_file, cond_file)
