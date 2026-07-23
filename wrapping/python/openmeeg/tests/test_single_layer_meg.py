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

import logging
import os

import numpy as np
import pytest
from numpy.testing import assert_array_less

import openmeeg as om

MAG_FACTOR = 1e-7  # mu_0 / (4 * pi)

_log = logging.getLogger(__name__)


def _log_blas_environment():
    """Log BLAS/LAPACK configuration, loaded libraries, and thread env vars.

    OpenMEEG's C++ core links against the *system* BLAS (found by CMake at
    build time), while numpy typically ships its own bundled BLAS. Both can
    end up loaded in the same process, each with its own thread pool -- a
    plausible source of CI flakiness (thread-pool oversubscription, and
    possibly different numerics from different BLAS implementations/versions
    for the same computation). This is emitted at INFO level so it shows up
    in pytest's "Captured log call" section when a test fails, without
    needing any special CI configuration.
    """
    try:
        config = np.show_config(mode="dicts")
        deps = config.get("Build Dependencies", {})
        blas = deps.get("blas", {})
        lapack = deps.get("lapack", {})
        _log.info(
            "numpy BLAS: %s %s | LAPACK: %s %s",
            blas.get("name"),
            blas.get("version"),
            lapack.get("name"),
            lapack.get("version"),
        )
    except Exception as exc:
        _log.info("Could not introspect numpy BLAS config: %s", exc)

    _log.info(
        "OMP_NUM_THREADS=%s OPENBLAS_NUM_THREADS=%s MKL_NUM_THREADS=%s",
        os.environ.get("OMP_NUM_THREADS"),
        os.environ.get("OPENBLAS_NUM_THREADS"),
        os.environ.get("MKL_NUM_THREADS"),
    )

    try:
        with open("/proc/self/maps") as fid:
            lines = fid.readlines()
    except OSError:
        _log.info("Loaded shared libraries: unavailable (non-Linux?)")
        return
    libs = set()
    for line in lines:
        if any(key in line.lower() for key in ("blas", "lapack", "mkl")):
            parts = line.split()
            if parts:
                libs.add(parts[-1])
    _log.info("Loaded BLAS/LAPACK/MKL shared libraries: %s", sorted(libs))


def _sym_to_dense(sym):
    """Reconstruct a dense array from an OpenMEEG SymMatrix.

    Vectorized (column-at-a-time) unpacking of the column-major
    upper-triangular packed storage used by ``SymMatrix::operator()``
    (``data()[i + j*(j+1)/2]`` for ``i <= j``), which is much faster than
    calling ``.value(i, j)`` for every entry.
    """
    n = sym.nlin()
    flat = sym.array_flat()
    A = np.empty((n, n))
    idx = 0
    for j in range(n):
        col = flat[idx : idx + j + 1]
        A[: j + 1, j] = col
        A[j, : j + 1] = col
        idx += j + 1
    return A


def _log_headmat_diagnostics(label, headmat, headmat_inv):
    """Log HeadMat conditioning and the accuracy of its computed inverse.

    Helps confirm/refute whether numerical ill-conditioning (as opposed to
    e.g. a BLAS-vendor-dependent inversion bug) explains observed RDM drift
    for a given configuration.
    """
    A = _sym_to_dense(headmat)
    Ainv = _sym_to_dense(headmat_inv)
    eigvals = np.linalg.eigvalsh(A)
    abs_eig = np.abs(eigvals)
    cond = abs_eig.max() / abs_eig.min()
    residual = np.linalg.norm(A @ Ainv - np.eye(A.shape[0])) / A.shape[0]
    _log.info(
        "[%s] HeadMat %dx%d: cond=%.3e, |eig| in [%.3e, %.3e], "
        "mean(|A @ Ainv - I|)=%.3e",
        label,
        A.shape[0],
        A.shape[0],
        cond,
        abs_eig.min(),
        abs_eig.max(),
        residual,
    )


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


def _om_meg_gain(meshes, conductivity, dipoles, squids_file, label):
    _log_blas_environment()
    geom = om.make_nested_geometry(meshes, conductivity=conductivity)
    integrator = om.Integrator(3, 0, 0.005)
    headmat = om.HeadMat(geom, integrator)
    hminv = headmat.inverse()
    _log_headmat_diagnostics(label, headmat, hminv)
    dip = om.Matrix(np.asfortranarray(dipoles))
    dsm = om.DipSourceMat(geom, dip, integrator, "Brain")
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
    "n_layers, conductivity, tol",
    [
        (1, (0.3,), 0.02),
        # Kept at main's original 0.02 tolerance (rather than loosened)
        # deliberately: this case has been observed to drift as high as
        # ~0.025 on CI, and checking HeadMat's condition number locally
        # shows it is *not* ill-conditioned (cond ~3e2, vs ~1e16 --
        # essentially singular, by design -- for the 1-layer case above,
        # which is not flaky). So conditioning alone does not explain the
        # drift. Leaving the tight tolerance in place lets CI reproduce the
        # failure with _log_blas_environment/_log_headmat_diagnostics
        # (emitted at INFO level, surfaced by pytest on failure) attached,
        # to gather real evidence from the flaky runner instead of guessing.
        (3, (0.3, 0.006, 0.3), 0.02),
    ],
)
def test_meg_sphere_vs_sarvas(n_layers, conductivity, tol, tmp_path, caplog):
    """MEG leadfield on a sphere matches the analytic Sarvas solution."""
    caplog.set_level(logging.INFO, logger=__name__)
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
    gain = _om_meg_gain(meshes, conductivity, dipoles, squids, f"{n_layers}-layer")

    rdm, mag = _rdm_mag(gain, ref)
    # A wrong-sign or missing secondary term (cf. issue #577) blows RDM up well
    # past this; a correct forward on this mesh gives RDM ~1e-3.
    assert_array_less(rdm, tol)
    assert_array_less(np.abs(mag - 1.0), tol)


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

    gain = _om_meg_gain(
        [_icosphere(2, radius)], (0.3,), dipoles, squids, "1-layer-radial"
    )
    radial_norm = np.linalg.norm(gain[:, 0])
    tangential_norm = np.linalg.norm(gain[:, 1])
    assert radial_norm < 0.02 * tangential_norm
