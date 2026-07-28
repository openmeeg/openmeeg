import numpy as np
import pytest
from numpy.testing import assert_allclose

from openmeeg._utils import read_geom_meshes, read_mesh
from openmeeg._viz import main, plot_files, read_files


@pytest.mark.parametrize(
    "geom_file",
    [
        "Head1/Head1.geom",  # 1.1, Interface: "..."
        "Head1/Head1.geom_v2",  # 1.1, Mesh: "..." + numbered interfaces
        "Head1/Head1_legacy.geom",  # 1.0, bare file names
    ],
)
def test_read_geom_meshes(data_path, geom_file):
    """All .geom flavors must yield the same meshes (gh-593)."""
    mesh_files = read_geom_meshes(data_path / geom_file)
    assert [mesh_file.name for mesh_file in mesh_files] == [
        "cortex.1.tri",
        "skull.1.tri",
        "scalp.1.tri",
    ]
    # named interfaces, and a different subject entirely
    mesh_files = read_geom_meshes(data_path / "mne_sample_ico3/mne_sample_ico3.geom")
    assert [mesh_file.name for mesh_file in mesh_files] == [
        "outer_skull.tri",
        "inner_skull.tri",
        "outer_skin.tri",
    ]


def test_read_geom_meshes_missing(data_path, tmp_path):
    geom_file = tmp_path / "bad.geom"
    geom_file.write_text('# Domain Description 1.1\n\nInterface: "nonexistent.tri"\n')
    with pytest.raises(FileNotFoundError, match="nonexistent.tri"):
        read_geom_meshes(geom_file)


def test_read_mesh(data_path):
    """read_mesh must agree with the C++ reader."""
    points, triangles = read_mesh(data_path / "Head1" / "cortex.1.tri")
    assert points.shape == (42, 3)
    assert triangles.shape == (80, 3)
    assert triangles.min() == 0
    assert triangles.max() == len(points) - 1
    with pytest.raises(ValueError, match="Unsupported mesh extension"):
        read_mesh(data_path / "Head1" / "Head1.cond")


def test_read_files(data_path):
    """Reading a full head model should give 3 surfaces and 2 point sets."""
    head1 = data_path / "Head1"
    meshes, point_sets = read_files(
        [head1 / "Head1.geom", head1 / "Head1.dip", head1 / "Head1.squids"]
    )
    assert len(meshes) == 3
    assert [len(points) for points, _ in meshes] == [42, 42, 42]

    assert len(point_sets) == 2
    dipoles, moments = point_sets[0]
    assert dipoles.shape == (6, 3)
    assert moments.shape == (6, 3)
    squids, orientations = point_sets[1]
    assert squids.shape == (162, 3)
    assert orientations.shape == (162, 3)
    # .squids files have a label column that np.loadtxt cannot handle
    assert_allclose(squids[0], [-0.63087720, 0.0, 1.0207812], rtol=1e-6)

    # sensors without orientations
    meshes, point_sets = read_files([head1 / "Head1.patches"])
    assert meshes == []
    positions, orientations = point_sets[0]
    assert positions.shape == (32, 3)
    assert orientations is None


def test_read_files_errors(data_path, tmp_path):
    with pytest.raises(ValueError, match="Unknown extension '.cond'"):
        read_files([data_path / "Head1" / "Head1.cond"])
    with pytest.raises(FileNotFoundError, match="does-not-exist.tri"):
        read_files([tmp_path / "does-not-exist.tri"])
    bad = tmp_path / "bad.txt"
    np.savetxt(bad, np.zeros((3, 2)))
    with pytest.raises(ValueError, match="at least 3 columns"):
        read_files([bad])


def test_om_viz_help(capsys):
    """om_viz -h must work without any plotting dependency installed."""
    with pytest.raises(SystemExit) as exc:
        main(["-h"])
    assert exc.value.code == 0
    assert "om_viz" in capsys.readouterr().out


def test_plot_files(data_path):
    """Smoke test the actual plotting."""
    pytest.importorskip("pyvistaqt")
    head1 = data_path / "Head1"
    plotter = plot_files(
        [head1 / "Head1.geom", head1 / "Head1.dip", head1 / "Head1.squids"],
        show=False,
        off_screen=True,
    )
    try:
        # 3 surfaces + (points, arrows) for the dipoles + same for the squids
        assert len(plotter.renderer.actors) == 7
        plotter.app.processEvents()
    finally:
        plotter.close()

    with pytest.raises(ValueError, match="No files to plot"):
        plot_files([], show=False, off_screen=True)
