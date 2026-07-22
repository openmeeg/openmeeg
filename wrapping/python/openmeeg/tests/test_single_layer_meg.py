"""Validate single- and multi-layer MEG forward against the analytic sphere.

Regression test for https://github.com/openmeeg/openmeeg/issues/577, which
reported that the 1-layer MEG leadfield disagreed with reference solutions
while the 3-layer model was fine (and noted that CTest had no 1-layer tests).

For a spherically symmetric conductor the external magnetic field is
independent of the conductivity profile (Sarvas, 1987), so a BEM sphere -- for
any number of layers -- must reproduce the analytic Sarvas field. We build
icospheres (regular, outward normals) rather than reusing the bundled Head0
meshes, whose coarse/inward-normal triangulation makes the 1-layer near-field
cancellation needlessly hard to resolve.
"""

import numpy as np
import pytest
from numpy.testing import assert_array_less

import openmeeg as om

MAG_FACTOR = 1e-7  # mu_0 / (4 * pi)


def _icosphere(n_sub, radius):
    """Subdivided icosahedron with outward-pointing triangle normals."""
    t = (1.0 + 5.0**0.5) / 2.0
    verts = [
        [-1, t, 0],
        [1, t, 0],
        [-1, -t, 0],
        [1, -t, 0],
        [0, -1, t],
        [0, 1, t],
        [0, -1, -t],
        [0, 1, -t],
        [t, 0, -1],
        [t, 0, 1],
        [-t, 0, -1],
        [-t, 0, 1],
    ]
    verts = [np.array(v, float) for v in verts]
    faces = [
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ]
    for _ in range(n_sub):
        mid = {}
        new_faces = []

        def _midpoint(a, b):
            key = (a, b) if a < b else (b, a)
            if key not in mid:
                verts.append((verts[a] + verts[b]) / 2.0)
                mid[key] = len(verts) - 1
            return mid[key]

        for a, b, c in faces:
            ab, bc, ca = _midpoint(a, b), _midpoint(b, c), _midpoint(c, a)
            new_faces += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        faces = new_faces
    verts = np.array(verts)
    verts *= radius / np.linalg.norm(verts, axis=1, keepdims=True)
    return verts, np.array(faces)


def _sarvas(rq, Q, r, ori):
    """Analytic magnetic field of a dipole Q at rq in a sphere at the origin.

    Returns the field projected onto ``ori`` at sensor location ``r``.
    """
    a = r - rq
    an = np.linalg.norm(a)
    rn = np.linalg.norm(r)
    F = an * (rn * an + rn**2 - np.dot(rq, r))
    grad_F = (an**2 / rn + np.dot(a, r) / an + 2 * an + 2 * rn) * r - (
        an + 2 * rn + np.dot(a, r) / an
    ) * rq
    Qxrq = np.cross(Q, rq)
    B = MAG_FACTOR / F**2 * (F * Qxrq - np.dot(Qxrq, r) * grad_F)
    return np.dot(B, ori)


def _om_meg_gain(meshes, conductivity, dipoles, squids_file):
    geom = om.make_nested_geometry(meshes, conductivity=conductivity)
    integrator = om.Integrator(3, 0, 0.005)
    hminv = om.HeadMat(geom, integrator).inverse()
    dip = om.Matrix(np.asfortranarray(dipoles))
    dsm = om.DipSourceMat(geom, dip, "Brain")
    sensors = om.Sensors(str(squids_file))
    h2mm = om.Head2MEGMat(geom, sensors)
    ds2mm = om.DipSource2MEGMat(dip, sensors)
    return om.GainMEG(hminv, dsm, h2mm, ds2mm).array()


def _rdm_mag(gain, ref):
    """Relative difference measure and magnitude ratio, per column (dipole)."""
    a = gain / np.linalg.norm(gain, axis=0, keepdims=True)
    b = ref / np.linalg.norm(ref, axis=0, keepdims=True)
    rdm = np.linalg.norm(a - b, axis=0)
    mag = np.linalg.norm(gain, axis=0) / np.linalg.norm(ref, axis=0)
    return rdm, mag


@pytest.mark.parametrize(
    "n_layers, conductivity",
    [
        (1, (0.3,)),
        (3, (0.3, 0.006, 0.3)),
    ],
)
def test_meg_sphere_vs_sarvas(n_layers, conductivity, tmp_path):
    """MEG leadfield on a sphere matches the analytic Sarvas solution."""
    radius = 0.1  # radius of the innermost (brain) sphere
    # Tangential dipoles along +z at increasing depth. Purely radial dipoles
    # produce no external field, so they are covered separately below.
    depths = np.array([0.5, 0.7, 0.8, 0.9]) * radius
    dipoles = np.zeros((len(depths), 6))
    dipoles[:, 2] = depths  # position on the z axis
    dipoles[:, 3] = 1.0  # moment along x (tangential)

    # Magnetometers on a ring outside the head, oriented radially.
    n_sensors = 32
    phi = np.linspace(0, 2 * np.pi, n_sensors, endpoint=False)
    theta = np.linspace(0.2, np.pi - 0.2, n_sensors)
    spos = (
        1.2
        * radius
        * np.c_[np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
    )
    sori = spos / np.linalg.norm(spos, axis=1, keepdims=True)
    squids = tmp_path / "sphere.squids"
    with open(squids, "w") as fid:
        for p, o in zip(spos, sori):
            fid.write("%g %g %g %g %g %g\n" % (*p, *o))

    # Analytic ground truth (sphere centered at the origin).
    ref = np.array(
        [
            [
                _sarvas(dipoles[d, :3], dipoles[d, 3:], spos[s], sori[s])
                for d in range(len(dipoles))
            ]
            for s in range(n_sensors)
        ]
    )

    # Three shells (brain/skull/scalp); only the inner one is used for 1 layer.
    meshes = [_icosphere(2, radius * f) for f in (1.0, 1.1, 1.2)][:n_layers]
    gain = _om_meg_gain(meshes, conductivity, dipoles, squids)

    rdm, mag = _rdm_mag(gain, ref)
    # A wrong-sign or missing secondary term (cf. issue #577) blows RDM up well
    # past this; a correct forward on this mesh gives RDM ~1e-3.
    print(rdm)
    assert_array_less(rdm, 0.02)
    print(np.abs(mag - 1.0))
    assert_array_less(np.abs(mag - 1.0), 0.02)


def test_meg_sphere_radial_dipole_is_silent(tmp_path):
    """A purely radial dipole produces (numerically) no external MEG field."""
    radius = 0.1
    # One radial (z) and one tangential (x) dipole at the same depth.
    dipoles = np.array(
        [
            [0.0, 0.0, 0.08, 0.0, 0.0, 1.0],  # radial -> silent
            [0.0, 0.0, 0.08, 1.0, 0.0, 0.0],  # tangential -> reference scale
        ]
    )
    n_sensors = 32
    phi = np.linspace(0, 2 * np.pi, n_sensors, endpoint=False)
    theta = np.linspace(0.2, np.pi - 0.2, n_sensors)
    spos = (
        1.2
        * radius
        * np.c_[np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
    )
    sori = spos / np.linalg.norm(spos, axis=1, keepdims=True)
    squids = tmp_path / "sphere.squids"
    with open(squids, "w") as fid:
        for p, o in zip(spos, sori):
            fid.write("%g %g %g %g %g %g\n" % (*p, *o))

    gain = _om_meg_gain([_icosphere(2, radius)], (0.3,), dipoles, squids)
    radial_norm = np.linalg.norm(gain[:, 0])
    tangential_norm = np.linalg.norm(gain[:, 1])
    assert radial_norm < 0.02 * tangential_norm
