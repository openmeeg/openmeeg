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
from os import path as op
import pytest
import openmeeg as om


@pytest.mark.skipif(os.getenv('OPENMEEG_BAD_MSVC') == '1',
                    reason="bug with msvc-based wrapping")
def test_python(data_path):
    # Load data
    subject = "Head1"
    cond_file = op.join(data_path, subject, subject + ".cond")
    geom_file = op.join(data_path, subject, subject + ".geom")
    source_mesh_file = op.join(data_path, subject, subject + ".tri")
    dipole_file = op.join(data_path, subject, subject + ".dip")
    squidsFile = op.join(data_path, subject, subject + ".squids")
    patches_file = op.join(data_path, subject, subject + ".patches")

    geom = om.Geometry(geom_file, cond_file)

    mesh = om.Mesh(source_mesh_file)

    dipoles = om.Matrix()
    dipoles.load(dipole_file)

    sensors = om.Sensors()
    sensors.load(squidsFile)

    patches = om.Sensors()
    patches.load(patches_file)

    # Compute forward problem (Build Gain Matrices)

    gauss_order = 3
    use_adaptive_integration = True
    dipole_in_cortex = True

    hm = om.HeadMat(geom)
    # hm.invert() # invert hm inplace (no copy)
    # hminv = hm
    hminv = hm.inverse()  # invert hm with a copy
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


    print("hm                  : %d x %d" % (hm.nlin(), hm.ncol()))
    print("hminv               : %d x %d" % (hminv.nlin(), hminv.ncol()))
    print("ssm                 : %d x %d" % (ssm.nlin(), ssm.ncol()))
    print("ss2mm               : %d x %d" % (ss2mm.nlin(), ss2mm.ncol()))
    print("dsm                 : %d x %d" % (ssm.nlin(), ssm.ncol()))
    print("ds2mm               : %d x %d" % (ss2mm.nlin(), ss2mm.ncol()))
    print("h2mm                : %d x %d" % (h2mm.nlin(), h2mm.ncol()))
    print("h2em                : %d x %d" % (h2mm.nlin(), h2mm.ncol()))
    print(
        "gain_meg_surf       : %d x %d"
        % (gain_meg_surf.nlin(), gain_meg_surf.ncol())
    )
    print(
        "gain_eeg_surf       : %d x %d"
        % (gain_eeg_surf.nlin(), gain_eeg_surf.ncol())
    )
    print(
        "gain_meg_dip        : %d x %d" % (gain_meg_dip.nlin(), gain_meg_dip.ncol())
    )
    print(
        "gain_adjoint_meg_dip: %d x %d"
        % (gain_adjoint_meg_dip.nlin(), gain_adjoint_meg_dip.ncol())
    )
    print(
        "gain_eeg_dip        : %d x %d" % (gain_eeg_dip.nlin(), gain_eeg_dip.ncol())
    )
    print(
        "gain_adjoint_eeg_dip: %d x %d"
        % (gain_adjoint_eeg_dip.nlin(), gain_adjoint_eeg_dip.ncol())
    )

    # Leadfield MEG in one line :

    gain_meg_surf_one_line = om.GainMEG(
        om.HeadMat(geom).inverse(),
        om.SurfSourceMat(geom, mesh),
        om.Head2MEGMat(geom, sensors),
        om.SurfSource2MEGMat(mesh, sensors),
    )

    print(
        "gain_meg_surf_one_line : %d x %d"
        % (gain_meg_surf_one_line.nlin(), gain_meg_surf_one_line.ncol())
    )

    ###############################################################################
    # Compute forward data =

    srcFile = op.join(data_path, subject, subject + ".srcdip")
    sources = om.Matrix(srcFile)

    noise_level = 0.0
    est_meg = om.Forward(gain_meg_dip, sources, noise_level)
    print("est_meg    : %d x %d" % (est_meg.nlin(), est_meg.ncol()))

    est_meg_adjoint = om.Forward(gain_adjoint_meg_dip, sources, noise_level)
    print(
        "est_meg_adjoint    : %d x %d"
        % (est_meg_adjoint.nlin(), est_meg_adjoint.ncol())
    )

    est_eeg = om.Forward(gain_eeg_dip, sources, noise_level)
    print("est_eeg    : %d x %d" % (est_eeg.nlin(), est_eeg.ncol()))

    est_eeg_adjoint = om.Forward(gain_adjoint_eeg_dip, sources, noise_level)
    print(
        "est_eeg_adjoint    : %d x %d"
        % (est_eeg_adjoint.nlin(), est_eeg_adjoint.ncol())
    )

    # Example of basic manipulations

    # TODO: the same with numpy
    v1 = om.Vertex(1.0, 0.0, 0.0, 0)
    v2 = om.Vertex(0.0, 1.0, 0.0, 1)
    v3 = om.Vertex(0.0, 0.0, 1.0, 2)
    # TODO: v4 = om.Vertex( [double] , int )

    # print(v1.norm()
    # print((v1 + v2).norm()

    normal = om.Vect3(1.0, 0.0, 0.0)
    t = om.Triangle(v1, v2, v3)

    hm_file = subject + ".hm"
    hm.save(hm_file)

    ssm_file = subject + ".ssm"
    ssm.save(ssm_file)

    m1 = om.SymMatrix()
    m1.load(hm_file)
    # print(m1(0, 0))
    # print(m1.nlin())
    # print(m1.ncol())

    m2 = om.Matrix()
    m2.load(ssm_file)
    # m2.setvalue(2,3,-0.2) # m2(2,3)=-0.2
    # print(m2(2,3))
    # print(m2(0, 0))
    # print(m2.nlin())
    # print(m2.ncol())

    # remove useless files
    os.remove(hm_file)
    os.remove(ssm_file)
