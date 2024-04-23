###########################################
# CAVEAT on using OpenMEEG from python !!
# Beware that a temporary object has its memory
# released. So do not work with data provided from
# an OpenMEEG temporary object.
# For example, having a symmetric matrix defined as :
# > M = om.SymMatrix(100)
#
# TODO:
# taking as an numpy array the sub-matrix of M might lead to corrupted memory:
# > mySubMat = om.asarray(M.submat(0,10,0,10))
# since submat returns a newly created object that is hold by nothing, thus
# destroyed afterward.
# ENDTODO
# Instead do keep an object pointing the newly created submatrix, and then
# access the numpy array form of it:
# > subM = M.submat(0,10,0,10)
# > mySubMat = om.asarray(subM)
###########################################
import os
import os.path as op
import shutil

import numpy as np
import pytest
from numpy.testing import assert_allclose

import openmeeg as om
import openmeeg._openmeeg_wrapper as _omc


def read_geom(geom_file):
    """ReadGeom : provides paths to meshes present in .geom file."""
    f = open(geom_file, "r")
    lines = f.readlines()
    mesh_files = []
    for line in lines:
        words = line.split()
        if len(words) > 1:
            if words[0] == "Interfaces":
                nb_mesh = int(words[1])
                print("Nb mesh files : %d" % nb_mesh)
                continue

        if len(words) == 1:
            mesh_file = words[0]
            if mesh_file.endswith(".tri"):
                mesh_files.append(mesh_file)
            if not os.path.exists(mesh_file):
                print("Could not find mesh : " + mesh_file)
            continue

        if (len(words) > 1) and words[0].startswith("Interface"):
            mesh_file = words[-1][1:-1]
            if mesh_file.endswith(".tri"):
                mesh_files.append(mesh_file)
            if not os.path.exists(mesh_file):
                print("Could not find mesh : " + mesh_file)
            continue

    for k, fname in enumerate(mesh_files):
        if not os.path.isabs(fname):
            mesh_files[k] = os.path.join(os.path.dirname(geom_file), fname)

    print("Found : %s" % mesh_files)
    return mesh_files


def read_tri(fname):
    """Read .tri file.

    Parameters
    ----------
    fname : str
        The file to read.

    Returns
    -------
    points : ndarray, shape (n_points, 3)
        The vertices
    normals : ndarray, shape (n_points, 3)
        The normals at the vertices
    faces : ndarray, shape (n_faces, 3)
        The faces
    """
    assert fname.endswith(".tri")
    fid = open(fname, "r")
    # read the number of vertices
    npoints = int(fid.readline().split()[1])

    points = []
    faces = []
    normals = []

    # fills the vertices arrays
    for _ in range(npoints):
        vals = fid.readline().split()
        points.append(vals[:3])
        normals.append(vals[3:])

    # Read the number of triangles
    n_faces = int(fid.readline().split()[1])
    # create the list of triangles
    for _ in range(n_faces):
        vals = fid.readline().split()
        faces.append(vals[:3])

    # Convert to numpy arrays
    points = np.array(points, np.float64)
    normals = np.array(normals, np.float64)
    faces = np.array(faces, np.int64)
    return points, normals, faces


@pytest.mark.parametrize(
    "subject",
    [
        "Head1",
        pytest.param("mne_sample_ico3", marks=pytest.mark.slow),
    ],
)
@pytest.mark.parametrize("load_from_numpy", [True, False])
def test_python(subject, data_path, load_from_numpy, tmp_path):
    # Load data
    cond_file = op.join(data_path, subject, subject + ".cond")
    geom_file = op.join(data_path, subject, subject + ".geom")
    source_mesh_file = op.join(data_path, subject, subject + ".tri")
    dipole_file = op.join(data_path, subject, subject + ".dip")
    squidsFile = op.join(data_path, subject, subject + ".squids")
    patches_file = op.join(data_path, subject, subject + ".patches")

    if load_from_numpy:
        mesh_files = read_geom(geom_file)
        meshes = []
        for fname in mesh_files:  # we assume the order is form inner to outer
            points, _, tris = read_tri(fname)
            meshes.append((points, tris))
        geom = om.make_nested_geometry(meshes, conductivity=(1, 0.0125, 1))
        dipoles = np.loadtxt(dipole_file)
        dipoles = om.Matrix(np.asfortranarray(dipoles))
        usecols = range(1, 4)
        if subject.startswith("mne_"):
            usecols = range(3)
        patches = np.loadtxt(patches_file, usecols=usecols)
        patches = om.Matrix(np.asfortranarray(patches))
        patches = om.Sensors(patches, geom)
    else:
        geom = om.read_geometry(geom_file, cond_file)
        dipoles = om.Matrix()
        dipoles.load(dipole_file)
        patches = om.Sensors()
        patches.load(patches_file)

    assert geom.is_nested()

    n_dipoles = dipoles.nlin()

    has_meg = False
    if op.exists(squidsFile):
        has_meg = True
        sensors = om.Sensors()
        sensors.load(squidsFile)
        n_meg_sensors = sensors.getNumberOfSensors()

    n_eeg_sensors = patches.getNumberOfSensors()

    # Compute forward problem (Build Gain Matrices)
    # gauss_order = 3  # XXX Integrator is now exposed, just need to use it...
    # use_adaptive_integration = True
    # dipole_in_cortex = True

    integrator = om.Integrator(3, 0, 0.005)
    print()
    print("*" * 80)
    print("HeadMat instantiation:")
    with om.use_log_level("debug"):
        hm = om.HeadMat(geom, integrator)
    hminv = hm.inverse()  # invert hm with a copy
    with om.use_log_level("debug"):
        hminv_inplace = om.HeadMat(geom, integrator)
    hminv_inplace.invert()  # invert hm inplace (no copy)
    assert_allclose(om.Matrix(hminv).array(), om.Matrix(hminv_inplace).array())

    dsm = om.DipSourceMat(geom, dipoles, "Brain")
    h2em = om.Head2EEGMat(geom, patches)

    gain_eeg_dip = om.GainEEG(hminv, dsm, h2em)
    gain_adjoint_eeg_dip = om.GainEEGadjoint(geom, dipoles, hm, h2em)
    n_hm_unknowns = geom.nb_parameters() - geom.nb_current_barrier_triangles()

    assert hm.nlin() == hm.ncol() == n_hm_unknowns
    assert hminv.nlin() == hminv.ncol() == n_hm_unknowns
    assert (dsm.nlin(), dsm.ncol()) == (n_hm_unknowns, n_dipoles)
    assert (gain_eeg_dip.nlin(), gain_eeg_dip.ncol()) == (n_eeg_sensors, n_dipoles)
    assert (gain_adjoint_eeg_dip.nlin(), gain_adjoint_eeg_dip.ncol()) == (
        n_eeg_sensors,
        n_dipoles,
    )

    if op.exists(source_mesh_file):
        mesh = _omc.Mesh(source_mesh_file)
        n_dipoles_surf = mesh.vertices().size()
        ssm = om.SurfSourceMat(geom, mesh, integrator)
        gain_eeg_surf = om.GainEEG(hminv, ssm, h2em)
        assert (ssm.nlin(), ssm.ncol()) == (n_hm_unknowns, n_dipoles_surf)
        assert (gain_eeg_surf.nlin(), gain_eeg_surf.ncol()) == (
            n_eeg_sensors,
            n_dipoles_surf,
        )
        assert (gain_eeg_surf.nlin(), gain_eeg_surf.ncol()) == (
            n_eeg_sensors,
            n_dipoles_surf,
        )

    if has_meg:
        ss2mm = om.SurfSource2MEGMat(mesh, sensors)
        ds2mm = om.DipSource2MEGMat(dipoles, sensors)
        # XXX test DipSourceMat with dipoles as array and make sure it breaks
        # if it's transposed
        h2mm = om.Head2MEGMat(geom, sensors)
        gain_meg_surf = om.GainMEG(hminv, ssm, h2mm, ss2mm)
        gain_meg_dip = om.GainMEG(hminv, dsm, h2mm, ds2mm)
        gain_adjoint_meg_dip = om.GainMEGadjoint(geom, dipoles, hm, h2mm, ds2mm)
        gain_adjoint_eeg_meg_dip = om.GainEEGMEGadjoint(
            geom, dipoles, hm, h2em, h2mm, ds2mm
        )

        assert gain_adjoint_eeg_meg_dip.nlin() == (n_meg_sensors + n_eeg_sensors)
        assert (ss2mm.nlin(), ss2mm.ncol()) == (n_meg_sensors, n_dipoles_surf)
        assert (ds2mm.nlin(), ds2mm.ncol()) == (n_meg_sensors, n_dipoles)
        assert (h2mm.nlin(), h2mm.ncol()) == (n_meg_sensors, n_hm_unknowns)
        assert (h2em.nlin(), h2em.ncol()) == (n_eeg_sensors, n_hm_unknowns)
        assert (gain_meg_dip.nlin(), gain_meg_dip.ncol()) == (n_meg_sensors, n_dipoles)
        assert (gain_meg_surf.nlin(), gain_meg_surf.ncol()) == (
            n_meg_sensors,
            n_dipoles_surf,
        )
        assert (gain_adjoint_meg_dip.nlin(), gain_adjoint_meg_dip.ncol()) == (
            n_meg_sensors,
            n_dipoles,
        )

        # Leadfield MEG in one line :
        gain_meg_surf_one_line = om.GainMEG(
            om.HeadMat(geom, integrator).inverse(),
            om.SurfSourceMat(geom, mesh, integrator),
            om.Head2MEGMat(geom, sensors),
            om.SurfSource2MEGMat(mesh, sensors),
        )

        assert_allclose(gain_meg_surf_one_line.array(), gain_meg_surf.array())

    if subject != "Head1":
        return  # only test the rest for Head1

    ###########################################################################
    # Compute forward data =

    src_fname = op.join(data_path, subject, subject + ".srcdip")
    sources = _omc.Matrix(src_fname)
    n_times = sources.ncol()

    noise_level = 0.0
    est_meg = om.Forward(gain_meg_dip, sources, noise_level)
    assert (est_meg.nlin(), est_meg.ncol()) == (n_meg_sensors, n_times)

    est_meg_adjoint = om.Forward(gain_adjoint_meg_dip, sources, noise_level)
    assert (est_meg_adjoint.nlin(), est_meg_adjoint.ncol()) == (n_meg_sensors, n_times)

    est_eeg = om.Forward(gain_eeg_dip, sources, noise_level)
    assert (est_eeg.nlin(), est_eeg.ncol()) == (n_eeg_sensors, n_times)

    est_eeg_adjoint = om.Forward(gain_adjoint_eeg_dip, sources, noise_level)
    assert (est_eeg_adjoint.nlin(), est_eeg_adjoint.ncol()) == (n_eeg_sensors, n_times)

    # Example of basic manipulations
    # TODO: the same with numpy
    v1 = _omc.Vertex(1.0, 0.0, 0.0, 0)
    v2 = _omc.Vertex(0.0, 1.0, 0.0, 1)
    _omc.Vertex(0.0, 0.0, 1.0, 2)
    # TODO: v4 = om.Vertex( [double] , int )

    # print(v1.norm()
    assert v1.norm() == np.linalg.norm(v1.array())
    assert (v1 + v2).norm() == np.linalg.norm(v1.array() + v2.array())

    normal = _omc.Vect3(1.0, 0.0, 0.0)
    assert normal.norm() == 1.0

    hm_fname = str(tmp_path / f"{subject}.hm")
    hm.save(hm_fname)

    ssm_fname = str(tmp_path / f"{subject}.ssm")
    ssm.save(ssm_fname)

    m1 = om.SymMatrix()
    m1.load(hm_fname)

    m2 = om.Matrix()
    m2.load(ssm_fname)

    assert_allclose(om.Matrix(hm).array(), om.Matrix(m1).array())
    assert_allclose(ssm.array(), m2.array())


def test_pointer_clearing(data_path, tmp_path):
    """Test that file pointers are cleared after data is loaded."""
    # copy files to tmp_path
    for fname in ("cortex.1.tri", "skull.1.tri", "scalp.1.tri"):
        shutil.copyfile(data_path / "Head1" / fname, tmp_path / fname)
    for ext in ("cond", "geom", "tri", "dip", "squids", "patches"):
        fname = f"Head1.{ext}"
        shutil.copyfile(data_path / "Head1" / fname, tmp_path / fname)
    del data_path

    # now load each one and delete it immediately
    tmp_geom_path = tmp_path / "Head1.geom"
    tmp_cond_path = tmp_path / "Head1.cond"
    assert tmp_geom_path.is_file()
    assert tmp_cond_path.is_file()
    geom = om.read_geometry(str(tmp_geom_path), str(tmp_cond_path))
    assert geom.is_nested()
    os.remove(tmp_geom_path)
    os.remove(tmp_cond_path)

    tmp_mesh_path = tmp_path / "Head1.tri"
    assert tmp_mesh_path.is_file()
    mesh = _omc.Mesh(str(tmp_mesh_path))
    assert mesh.vertices().size() == 42
    os.remove(tmp_mesh_path)

    tmp_dipole_path = tmp_path / "Head1.dip"
    assert tmp_dipole_path.is_file()
    dipoles = om.Matrix()
    dipoles.load(str(tmp_dipole_path))
    assert dipoles.nlin() == 6
    os.remove(tmp_dipole_path)

    tmp_squids_path = tmp_path / "Head1.squids"
    assert tmp_squids_path.is_file()
    sensors = om.Sensors()
    sensors.load(str(tmp_squids_path))
    assert sensors.getNumberOfSensors() == 162

    tmp_patches_path = tmp_path / "Head1.patches"
    assert tmp_patches_path.is_file()
    patches = om.Sensors()
    patches.load(str(tmp_patches_path))
    assert patches.getNumberOfSensors() == 32
