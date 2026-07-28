"""Simple viewer for OpenMEEG head model files.

Originally written by Alexandre Gramfort <alexandre.gramfort@inria.fr>.
"""

import argparse
from pathlib import Path

import numpy as np

from ._openmeeg_wrapper import Sensors
from ._utils import MESH_EXTENSIONS, read_geom_meshes, read_mesh

# Sensor files are read through OpenMEEG's own reader, which knows how to deal
# with the optional label / weight / radius columns.
SENSOR_EXTENSIONS = (".electrodes", ".patches", ".squids")
# Plain (n_points, 3) or (n_points, 6) text files: positions [+ orientations]
POINT_EXTENSIONS = (".dip", ".src", ".txt")
ALL_EXTENSIONS = tuple(
    sorted(MESH_EXTENSIONS + SENSOR_EXTENSIONS + POINT_EXTENSIONS + (".geom",))
)

# light pink, khaki, light blue, grey
_MESH_COLORS = (
    (0.95, 0.83, 0.83),
    (0.91, 0.89, 0.67),
    (0.67, 0.89, 0.91),
    (0.68, 0.68, 0.68),
)
_POINT_COLOR = (1.0, 0.0, 0.0)
_ARROW_COLOR = (1.0, 0.0, 1.0)


def read_files(fnames):
    """Read OpenMEEG files, guessing the content from the extension.

    Parameters
    ----------
    fnames : list of path-like
        The files to read. ``.geom`` files are expanded into the mesh files
        that they reference.

    Returns
    -------
    meshes : list of tuple
        One ``(points, triangles)`` tuple per surface.
    point_sets : list of tuple
        One ``(positions, orientations)`` tuple per sensor or dipole file,
        where ``orientations`` is None if the file does not contain any.
    """
    meshes = list()
    point_sets = list()
    for fname in fnames:
        fname = Path(fname)
        if not fname.is_file():
            raise FileNotFoundError(f"File not found: {fname}")
        suffix = fname.suffix
        if suffix == ".geom":
            meshes.extend(read_mesh(mesh_file) for mesh_file in read_geom_meshes(fname))
        elif suffix in MESH_EXTENSIONS:
            meshes.append(read_mesh(fname))
        elif suffix in SENSOR_EXTENSIONS:
            sensors = Sensors()
            sensors.load(str(fname))
            orientations = None
            if sensors.hasOrientations():
                orientations = sensors.getOrientations().array()
            point_sets.append((sensors.getPositions().array(), orientations))
        elif suffix in POINT_EXTENSIONS:
            data = np.loadtxt(fname, ndmin=2)
            if data.shape[1] < 3:
                raise ValueError(
                    f"Expected at least 3 columns in {fname}, got {data.shape[1]}"
                )
            orientations = data[:, 3:6] if data.shape[1] >= 6 else None
            point_sets.append((data[:, :3], orientations))
        else:
            raise ValueError(
                f"Unknown extension {repr(suffix)} for {fname}, must be one of "
                f"{list(ALL_EXTENSIONS)}"
            )
    return meshes, point_sets


def _to_poly_data(points, triangles):
    import pyvista as pv

    faces = np.c_[np.full(len(triangles), 3, np.int64), triangles]
    return pv.PolyData(points, faces.ravel())


def _scale_factor(meshes, point_sets):
    """Get a glyph size suitable for the extent of the data."""
    all_points = [mesh[0] for mesh in meshes] + [points[0] for points in point_sets]
    all_points = np.concatenate(all_points, axis=0)
    extent = all_points.max(axis=0) - all_points.min(axis=0)
    return 0.05 * np.linalg.norm(extent)


def plot_files(fnames, *, opacity=0.4, show=True, off_screen=None):
    """Plot OpenMEEG files in a :class:`pyvistaqt.BackgroundPlotter`.

    Parameters
    ----------
    fnames : list of path-like
        The files to plot, see :func:`read_files`.
    opacity : float
        The opacity to use for the surfaces.
    show : bool
        Whether to show the plotter window.
    off_screen : bool | None
        Whether to render off screen (useful for testing).

    Returns
    -------
    plotter : instance of pyvistaqt.BackgroundPlotter
        The plotter.
    """
    import pyvista as pv
    from pyvistaqt import BackgroundPlotter

    meshes, point_sets = read_files(fnames)
    if not meshes and not point_sets:
        raise ValueError("No files to plot")
    factor = _scale_factor(meshes, point_sets)

    plotter = BackgroundPlotter(
        show=show, off_screen=off_screen, title="OpenMEEG", menu_bar=False, editor=False
    )
    for ii, (points, triangles) in enumerate(meshes):
        plotter.add_mesh(
            _to_poly_data(points, triangles),
            color=_MESH_COLORS[ii % len(_MESH_COLORS)],
            opacity=opacity,
        )
    for positions, orientations in point_sets:
        point_cloud = pv.PolyData(positions)
        plotter.add_mesh(
            point_cloud,
            color=_POINT_COLOR,
            point_size=6,
            render_points_as_spheres=True,
        )
        if orientations is not None:
            point_cloud["vectors"] = orientations
            arrows = point_cloud.glyph(orient="vectors", scale=False, factor=factor)
            plotter.add_mesh(arrows, color=_ARROW_COLOR)
    plotter.add_axes()
    plotter.reset_camera()
    return plotter


def main(argv=None):
    """Run the ``om_viz`` command-line tool."""
    parser = argparse.ArgumentParser(
        prog="om_viz",
        description=(
            "View OpenMEEG head model files. The content of each file is "
            "guessed from its extension, and a .geom file additionally pulls "
            "in the meshes that it references."
        ),
        epilog=(
            "example: om_viz model.geom mesh.tri dipoles.dip sensors.squids\n"
            f"known extensions: {' '.join(ALL_EXTENSIONS)}"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("fnames", nargs="+", metavar="FILE", help="file(s) to view")
    parser.add_argument(
        "--opacity",
        type=float,
        default=0.4,
        help="opacity of the surfaces (default: %(default)s)",
    )
    args = parser.parse_args(argv)
    plotter = plot_files(args.fnames, opacity=args.opacity)
    # Qt bindings disagree about exec vs exec_
    exec_ = getattr(plotter.app, "exec", None) or plotter.app.exec_
    exec_()


if __name__ == "__main__":
    main()
