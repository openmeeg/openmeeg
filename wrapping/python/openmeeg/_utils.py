import re
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from ._openmeeg_wrapper import DEBUG, ERROR, INFORMATION, WARNING, Logger, Mesh

_warn_map = dict(
    error=ERROR,
    warning=WARNING,
    info=INFORMATION,
    debug=DEBUG,
)


def set_log_level(level):
    """Set the logging level.

    Parameters
    ----------
    level : str
        Can be ``'error'``, ``'warning'``, ``'info'``, or ``'debug'``.
    """
    if not isinstance(level, str) or level not in _warn_map:
        raise ValueError(
            f"Unknown level {repr(level)}, must be one of {list(_warn_map)}"
        )
    Logger.logger().set_info_level(_warn_map[level])


def get_log_level():
    """Get the current logging level.

    Returns
    -------
    level : str
        See :func:`set_log_level`.
    """
    rev_map = {val: key for key, val in _warn_map.items()}
    return rev_map[Logger.logger().get_info_level()]


@contextmanager
def use_log_level(level):
    """Context manager for logging level.

    Parameters
    ----------
    level : str
        See :func:`set_log_level`.
    """
    old_level = get_log_level()
    set_log_level(level)
    try:
        yield
    finally:
        set_log_level(old_level)


# Extensions understood by OpenMEEG's C++ mesh readers (see OpenMEEG/include/MeshIOs)
MESH_EXTENSIONS = (".bnd", ".gii", ".mesh", ".off", ".tri", ".vtk")

_QUOTED = re.compile(r'"([^"]*)"')


def mesh_vertices_and_triangles(mesh):
    """Extract the vertices and triangles of a Mesh as NumPy arrays.

    Parameters
    ----------
    mesh : instance of Mesh
        The mesh.

    Returns
    -------
    vertices : ndarray, shape (n_vertices, 3)
        The vertices.
    triangles : ndarray, shape (n_triangles, 3)
        The triangles, as indices into ``vertices``.
    """
    mesh_vertices = mesh.geometry().vertices()
    vertices = np.array([vertex.array() for vertex in mesh_vertices], np.float64)
    mesh_triangles = mesh.triangles()
    triangles = np.array(
        [mesh.triangle(triangle).array() for triangle in mesh_triangles],
        dtype=np.int64,
    )
    return vertices, triangles


def read_tri(fname):
    """Read a ``.tri`` file.

    Parameters
    ----------
    fname : path-like
        The file to read.

    Returns
    -------
    points : ndarray, shape (n_points, 3)
        The vertices.
    normals : ndarray, shape (n_points, 3)
        The normals at the vertices.
    triangles : ndarray, shape (n_triangles, 3)
        The triangles, as indices into ``points``.
    """
    fname = Path(fname)
    assert fname.suffix == ".tri", fname
    with open(fname, "r") as fid:
        # read the number of vertices
        n_points = int(fid.readline().split()[1])

        points = []
        normals = []
        triangles = []

        # fill the vertices arrays
        for _ in range(n_points):
            vals = fid.readline().split()
            points.append(vals[:3])
            normals.append(vals[3:])

        # read the number of triangles, then the triangles themselves
        n_triangles = int(fid.readline().split()[1])
        for _ in range(n_triangles):
            triangles.append(fid.readline().split()[:3])

    points = np.array(points, np.float64)
    normals = np.array(normals, np.float64)
    triangles = np.array(triangles, np.int64)
    return points, normals, triangles


def read_mesh(fname):
    """Read a mesh file.

    Parameters
    ----------
    fname : path-like
        The file to read. Any format supported by OpenMEEG can be used, see
        ``MESH_EXTENSIONS``.

    Returns
    -------
    points : ndarray, shape (n_points, 3)
        The vertices.
    triangles : ndarray, shape (n_triangles, 3)
        The triangles, as indices into ``points``.
    """
    fname = Path(fname)
    if fname.suffix == ".tri":  # pure-Python reader, no C++ round-trip
        points, _, triangles = read_tri(fname)
        return points, triangles
    if fname.suffix not in MESH_EXTENSIONS:
        raise ValueError(
            f"Unsupported mesh extension {repr(fname.suffix)}, must be one of "
            f"{list(MESH_EXTENSIONS)}: {fname}"
        )
    mesh = Mesh()
    mesh.load(str(fname), False)  # verbose=False
    return mesh_vertices_and_triangles(mesh)


def read_geom_meshes(fname):
    """Get the paths of the mesh files referenced by a ``.geom`` file.

    Parameters
    ----------
    fname : path-like
        The ``.geom`` file to parse. Version 1.0 (legacy) and 1.1 files, with
        either ``Interface:`` or ``Mesh:`` declarations, are supported.

    Returns
    -------
    mesh_files : list of Path
        The mesh files, in the order in which they are declared, resolved
        relative to the directory containing ``fname``.
    """
    fname = Path(fname)
    mesh_files = []
    for line in fname.read_text().splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        # Quoted names ('Interface: "cortex.1.tri"', 'Mesh Cortex: "..."'), else
        # a bare file name on its own line (version 1.0 files).
        candidates = _QUOTED.findall(line) or ([line] if len(line.split()) == 1 else [])
        mesh_files.extend(
            candidate
            for candidate in candidates
            if Path(candidate).suffix in MESH_EXTENSIONS
        )
    mesh_files = [fname.parent / mesh_file for mesh_file in mesh_files]
    missing = [str(mesh_file) for mesh_file in mesh_files if not mesh_file.is_file()]
    if missing:
        raise FileNotFoundError(
            f"Mesh file(s) referenced by {fname} not found: {missing}"
        )
    return mesh_files
