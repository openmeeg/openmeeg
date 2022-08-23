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
import os.path as op
import numpy as np

from numpy.testing import assert_allclose
import openmeeg as om


def test_python(data_path, tmp_path):
    # Load data
    subject = "Head1"
    cond_file = op.join(data_path, subject, subject + ".cond")
    geom_file = op.join(data_path, subject, subject + ".geom")
    source_mesh_file = op.join(data_path, subject, subject + ".tri")
    dipole_file = op.join(data_path, subject, subject + ".dip")
    squidsFile = op.join(data_path, subject, subject + ".squids")
    patches_file = op.join(data_path, subject, subject + ".patches")

    geom = om.Geometry(geom_file, cond_file)
    assert geom.is_nested()

    mesh = om.Mesh(source_mesh_file)
    n_dipoles_surf = mesh.vertices().size()

    dipoles = om.Matrix()
    dipoles.load(dipole_file)
    n_dipoles = dipoles.nlin()

    sensors = om.Sensors()
    sensors.load(squidsFile)
    n_meg_sensors = sensors.getNumberOfSensors()

    patches = om.Sensors()
    patches.load(patches_file)
    n_eeg_sensors = patches.getNumberOfSensors()

    # Compute forward problem (Build Gain Matrices)
    # gauss_order = 3  # XXX cannot get Integrator exposed
    # use_adaptive_integration = True
    # dipole_in_cortex = True

    hm = om.HeadMat(geom)
    hminv = hm.inverse()  # invert hm with a copy
    hminv_inplace = om.HeadMat(geom)
    hminv_inplace.invert()  # invert hm inplace (no copy)
    assert_allclose(om.Matrix(hminv).array(), om.Matrix(hminv_inplace).array())

    ssm = om.SurfSourceMat(geom, mesh)
    ss2mm = om.SurfSource2MEGMat(mesh, sensors)
    dsm = om.DipSourceMat(geom, dipoles)
    ds2mm = om.DipSource2MEGMat(dipoles, sensors)
    h2mm = om.Head2MEGMat(geom, sensors)
    h2em = om.Head2EEGMat(geom, patches)
    gain_meg_surf = om.GainMEG(hminv, ssm, h2mm, ss2mm)
    gain_eeg_surf = om.GainEEG(hminv, ssm, h2em)
    gain_meg_dip = om.GainMEG(hminv, dsm, h2mm, ds2mm)
    gain_adjoint_meg_dip = om.GainMEGadjoint(geom, dipoles, hm, h2mm, ds2mm)
    gain_eeg_dip = om.GainEEG(hminv, dsm, h2em)
    gain_adjoint_eeg_dip = om.GainEEGadjoint(geom, dipoles, hm, h2em)
    gain_adjoint_eeg_meg_dip = om.GainEEGMEGadjoint(
        geom, dipoles, hm, h2em, h2mm, ds2mm
    )

    n_hm_unknowns = geom.nb_parameters() - geom.nb_current_barrier_triangles()
    assert gain_adjoint_eeg_meg_dip.nlin() == (n_meg_sensors + n_eeg_sensors)
    assert hm.nlin() == hm.ncol() == n_hm_unknowns
    assert hminv.nlin() == hminv.ncol() == n_hm_unknowns
    assert (ssm.nlin(), ssm.ncol()) == (n_hm_unknowns, n_dipoles_surf)
    assert (dsm.nlin(), dsm.ncol()) == (n_hm_unknowns, n_dipoles)
    assert (ss2mm.nlin(), ss2mm.ncol()) == (n_meg_sensors, n_dipoles_surf)
    assert (ds2mm.nlin(), ds2mm.ncol()) == (n_meg_sensors, n_dipoles)
    assert (h2mm.nlin(), h2mm.ncol()) == (n_meg_sensors, n_hm_unknowns)
    assert (h2em.nlin(), h2em.ncol()) == (n_eeg_sensors, n_hm_unknowns)
    assert (gain_eeg_dip.nlin(), gain_eeg_dip.ncol()) == (n_eeg_sensors, n_dipoles)
    assert (gain_eeg_surf.nlin(), gain_eeg_surf.ncol()) == (
        n_eeg_sensors,
        n_dipoles_surf,
    )
    assert (gain_meg_dip.nlin(), gain_meg_dip.ncol()) == (n_meg_sensors, n_dipoles)
    assert (gain_meg_surf.nlin(), gain_meg_surf.ncol()) == (
        n_meg_sensors,
        n_dipoles_surf,
    )
    assert (gain_adjoint_meg_dip.nlin(), gain_adjoint_meg_dip.ncol()) == (
        n_meg_sensors,
        n_dipoles,
    )
    assert (gain_adjoint_eeg_dip.nlin(), gain_adjoint_eeg_dip.ncol()) == (
        n_eeg_sensors,
        n_dipoles,
    )

    # Leadfield MEG in one line :
    gain_meg_surf_one_line = om.GainMEG(
        om.HeadMat(geom).inverse(),
        om.SurfSourceMat(geom, mesh),
        om.Head2MEGMat(geom, sensors),
        om.SurfSource2MEGMat(mesh, sensors),
    )

    assert_allclose(gain_meg_surf_one_line.array(), gain_meg_surf.array())

    ###############################################################################
    # Compute forward data =

    src_fname = op.join(data_path, subject, subject + ".srcdip")
    sources = om.Matrix(src_fname)
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
    v1 = om.Vertex(1.0, 0.0, 0.0, 0)
    v2 = om.Vertex(0.0, 1.0, 0.0, 1)
    om.Vertex(0.0, 0.0, 1.0, 2)
    # TODO: v4 = om.Vertex( [double] , int )

    # print(v1.norm()
    assert v1.norm() == np.linalg.norm(v1.array())
    assert (v1 + v2).norm() == np.linalg.norm(v1.array() + v2.array())

    normal = om.Vect3(1.0, 0.0, 0.0)
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
