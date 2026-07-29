#!/usr/bin/env python
"""Dump every intermediate of the 3-layer Sarvas MEG test for CI forensics.

Background: ``test_meg_sphere_vs_sarvas[3-conductivity1]`` intermittently fails
on the ``ubuntu-22.04 / OpenBLAS / pypy-3.11`` job with a *bit-identical*
``rdm = [0.024518, 0.020993, 0.018968, 0.014939]`` while the very same job
passes on other runs -- i.e. two deterministic regimes tied to the runner, not
float jitter.  Local work (Ubuntu 26.04 host and an ubuntu-22.04 container with
gcc 11.4 / OpenBLAS 0.3.20 / LAPACKE 3.10, CPython and PyPy 3.11.15) always
reproduces the *passing* regime bit-for-bit, and ruled out interpreter, OpenMP
and OpenBLAS thread counts, mesh orientation, and rounding in general.

So this script exists to answer, from the runner itself: which operator is the
first to differ, and what is different about that machine.  Run it on as many
runners as possible and diff the resulting ``.npz``/``.json`` pairs.

Usage:  python debug_sarvas_dump.py <output_dir> [label]

Never raises: a dump failure must not mask the thing we are trying to observe.
"""

from __future__ import annotations

import ctypes
import hashlib
import json
import os
import platform
import subprocess
import sys
import traceback

import numpy as np

MAG_FACTOR = 1e-7  # mu_0 / (4 * pi)

# The values seen on the failing runners, for quick eyeballing of the log.
CI_FAILING_RDM = [0.024518380512803167, 0.020992679522427646, 0.018968, 0.014939]


# --------------------------------------------------------------------------
# geometry / analytic reference (kept byte-identical to the test on purpose)
# --------------------------------------------------------------------------


def icosphere(n_sub, radius):
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


def sarvas(rq, Q, r, ori):
    """Analytic field of dipole Q at rq, projected on ``ori`` at ``r``."""
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


def rdm_mag(gain, ref):
    """Relative difference measure and magnitude ratio, per dipole."""
    a = gain / np.linalg.norm(gain, axis=0, keepdims=True)
    b = ref / np.linalg.norm(ref, axis=0, keepdims=True)
    rdm = np.linalg.norm(a - b, axis=0)
    mag = np.linalg.norm(gain, axis=0) / np.linalg.norm(ref, axis=0)
    return rdm, mag


def unpack_sym(sym):
    """Dense array from OpenMEEG's column-major packed upper-triangular store."""
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


def digest(a):
    """SHA-256 of an array's float64 bytes, for cheap cross-runner comparison."""
    a = np.ascontiguousarray(np.asarray(a, np.float64))
    return hashlib.sha256(a.tobytes()).hexdigest()


# --------------------------------------------------------------------------
# environment probes
# --------------------------------------------------------------------------


def _run(cmd):
    """Best-effort shell capture; diagnostics must never raise."""
    try:
        return subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=60
        ).stdout.strip()
    except Exception as exc:  # pragma: no cover - diagnostics only
        return f"<failed: {exc}>"


def loaded_blas_libs():
    """Paths of every BLAS/LAPACK/MKL shared object mapped into this process."""
    libs = set()
    try:
        with open("/proc/self/maps") as fid:
            for line in fid:
                low = line.lower()
                if any(k in low for k in ("blas", "lapack", "mkl")):
                    libs.add(line.split()[-1])
    except OSError:
        pass
    return sorted(libs)


def openblas_details(libs):
    """Ask each loaded OpenBLAS which kernel it picked -- the key CPU signal."""
    out = {}
    for lib in libs:
        if "openblas" not in os.path.basename(lib).lower():
            continue
        info = {}
        try:
            dll = ctypes.CDLL(lib)
        except OSError as exc:
            out[lib] = {"error": str(exc)}
            continue

        # Plain OpenBLAS exports openblas_get_config; the scipy-openblas wheels
        # rename to scipy_openblas_get_config64_ (prefix *and* ILP64 suffix).
        def _lookup(sym, restype, cast):
            for prefix in ("", "scipy_"):
                for suffix in ("", "64_", "_64_"):
                    try:
                        fn = getattr(dll, f"{prefix}{sym}{suffix}")
                    except AttributeError:
                        continue
                    fn.restype = restype
                    try:
                        info[sym] = cast(fn())
                    except Exception as exc:
                        info[sym] = f"<failed: {exc}>"
                    return
            info[sym] = "<symbol not found>"

        for sym in ("openblas_get_config", "openblas_get_corename"):
            _lookup(sym, ctypes.c_char_p, lambda v: v.decode(errors="replace"))
        for sym in ("openblas_get_num_threads", "openblas_get_num_procs"):
            _lookup(sym, ctypes.c_int, int)
        out[lib] = info
    return out


def binary_hashes(ldd_output, ext_paths):
    """SHA-256 of the OpenMEEG binaries actually loaded.

    HeadMat has now come back bit-identical on 20 runners across three
    different OpenBLAS 0.3.20 kernels (Zen, SkylakeX, Prescott), and the
    inverse is clean at every thread count.  Since the assembly uses no BLAS
    at all, a run that produces a different HeadMat would have to be running a
    *different binary* -- which is plausible because the workflow restores a
    ccache keyed only on job name and OS, shared across the whole matrix.  A
    stale or mismatched cache entry yielding a subtly different libOpenMEEG
    would be per-run rather than per-runner, and would explain values that are
    bit-identical across unrelated PRs yet come and go over time.
    """
    out = {}
    paths = list(ext_paths)
    for line in (ldd_output or "").splitlines():
        # "libOpenMEEG.so.1 => /path/to/libOpenMEEG.so.1.1.0 (0x...)"
        if "=>" in line and ("OpenMEEG" in line or "openblas" in line.lower()):
            target = line.split("=>", 1)[1].strip().split(" (")[0].strip()
            if target and os.path.isfile(target):
                paths.append(target)
    for path in paths:
        try:
            with open(path, "rb") as fid:
                out[os.path.basename(path)] = {
                    "sha256": hashlib.sha256(fid.read()).hexdigest(),
                    "size": os.path.getsize(path),
                    "path": path,
                }
        except OSError as exc:
            out[os.path.basename(path)] = {"error": str(exc)}
    return out


def openmeeg_openblas():
    """Ctypes handle for the OpenBLAS OpenMEEG links (not numpy's bundled one)."""
    for lib in loaded_blas_libs():
        base = os.path.basename(lib)
        if "openblas" in base.lower() and "scipy_openblas" not in base:
            try:
                return ctypes.CDLL(lib), lib
            except OSError:
                pass
    return None, None


def inverse_thread_sweep(headmat, A, cond):
    """Re-invert HeadMat at several OpenBLAS thread counts.

    OpenBLAS#3044 was a *wrong result* (failed convergence) specific to both an
    architecture and a thread count.  The analogue here is DSPTRF/DSPTRI
    returning Info==0 with a silently wrong factorization.  Because RDM depends
    on HeadMatInv linearly with gain 1.0, an inverse wrong by ~4x is exactly
    what would take max RDM from 0.0061 to 0.0245 -- so this is the one story
    that survives everything else we have excluded.

    Sweeping the thread count inside every job probes that failure mode
    deterministically, instead of waiting for a rare scheduling coincidence.
    """
    dll, path = openmeeg_openblas()
    out = {"library": path}
    if dll is None:
        out["error"] = "OpenMEEG's libopenblas not found among loaded libraries"
        return out
    try:
        setter = dll.openblas_set_num_threads
        getter = dll.openblas_get_num_threads
    except AttributeError as exc:
        out["error"] = f"no openblas_set_num_threads: {exc}"
        return out
    setter.argtypes = [ctypes.c_int]
    getter.restype = ctypes.c_int
    original = int(getter())
    n = A.shape[0]
    eye = np.eye(n)
    try:
        for nthread in (1, 2, 4, 8):
            setter(nthread)
            inv = unpack_sym(headmat.inverse())
            out[str(nthread)] = {
                "requested": nthread,
                "effective": int(getter()),
                "residual": float(np.linalg.norm(A @ inv - eye) / n),
                "digest": digest(inv),
            }
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        setter(original)
    resids = [
        v["residual"] for v in out.values() if isinstance(v, dict) and "residual" in v
    ]
    if resids:
        out["worst_residual"] = max(resids)
        out["residual_spread"] = max(resids) / (min(resids) or 1.0)
        # Do NOT flag purely on a large residual: the 1-layer HeadMat is
        # singular by construction (deflated current barrier, cond ~1e15), so
        # its inverse is meaningless and its residual is ~0.1 on every machine.
        # A real DSPTRF/DSPTRI failure shows up either as a bad residual on a
        # well-conditioned matrix, or as a residual that depends on the thread
        # count -- the OpenBLAS#3044 signature.
        out["any_broken"] = bool(
            (cond < 1e8 and max(resids) > 1e-10) or out["residual_spread"] > 1e3
        )
    return out


def environment():
    """Collect everything identifying about this runner and its BLAS."""
    # Import openmeeg FIRST.  Otherwise the system libopenblas that OpenMEEG
    # actually links is not yet mapped, and we would only ever see numpy's
    # bundled libscipy_openblas -- which is the wrong library for the question
    # we are asking (numpy's is a different OpenBLAS version and may dispatch
    # to a different kernel on the same CPU).
    try:
        import openmeeg  # noqa: F401
    except Exception:
        pass

    env = {
        "python_version": sys.version,
        "python_implementation": sys.implementation.name,
        "numpy_version": np.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "libc": "-".join(x for x in platform.libc_ver() if x),
        "env_vars": {
            k: v
            for k, v in sorted(os.environ.items())
            if "THREAD" in k or "OMP" in k or "BLAS" in k or "PYTHON" in k
        },
        "loaded_blas_libs": loaded_blas_libs(),
        "lscpu": _run("lscpu"),
        "cpu_model": _run("grep -m1 'model name' /proc/cpuinfo"),
        "cpu_flags": _run("grep -m1 '^flags' /proc/cpuinfo"),
        "gcc_version": _run("gcc --version | head -1"),
        "cpu_count": os.cpu_count(),
        "nproc": _run("nproc"),
        "ldd_openmeeg": "",
    }
    # OpenMEEG's assembly is OpenMP-parallel (libgomp), which is a different
    # thread pool from OpenBLAS's -- record it separately.
    try:
        # This probe runs before openmeeg is imported, so libgomp may not be
        # mapped yet -- fall back to the soname, which loads it either way.
        gomp = [
            line.split()[-1] for line in open("/proc/self/maps") if "libgomp" in line
        ]
        lib = ctypes.CDLL(gomp[0] if gomp else "libgomp.so.1")
        lib.omp_get_max_threads.restype = ctypes.c_int
        env["omp_get_max_threads"] = int(lib.omp_get_max_threads())
    except Exception as exc:
        env["omp_get_max_threads"] = f"<unavailable: {exc}>"
    env["openblas"] = openblas_details(env["loaded_blas_libs"])
    try:
        import openmeeg

        env["openmeeg_version"] = getattr(openmeeg, "__version__", "?")
        env["openmeeg_file"] = openmeeg.__file__
        so = [
            os.path.join(os.path.dirname(openmeeg.__file__), f)
            for f in os.listdir(os.path.dirname(openmeeg.__file__))
            if f.startswith("_openmeeg") and (".so" in f or ".pyd" in f)
        ]
        env["openmeeg_ext"] = so
        if so:
            env["ldd_openmeeg"] = _run(f"ldd {so[0]}")
        env["binary_hashes"] = binary_hashes(env["ldd_openmeeg"], so)
        env["ccache_version"] = _run("ccache --version 2>/dev/null | head -1")
        env["ccache_stats"] = _run("ccache -s 2>/dev/null")
    except Exception as exc:
        env["openmeeg_error"] = str(exc)
    try:
        import threadpoolctl

        env["threadpoolctl"] = threadpoolctl.threadpool_info()
    except Exception as exc:
        env["threadpoolctl"] = f"<unavailable: {exc}>"
    try:
        env["numpy_config"] = np.show_config(mode="dicts")
    except Exception as exc:
        env["numpy_config"] = f"<unavailable: {exc}>"
    return env


# --------------------------------------------------------------------------
# the computation itself
# --------------------------------------------------------------------------


def compute(n_layers, out, tag):
    """Run one forward and record every intermediate."""
    import openmeeg as om

    radius = 0.1
    depths = np.array([0.5, 0.7, 0.8, 0.9]) * radius
    dipoles = np.zeros((len(depths), 6))
    dipoles[:, 2] = depths
    dipoles[:, 3] = 1.0

    n_sensors = 32
    phi = np.linspace(0, 2 * np.pi, n_sensors, endpoint=False)
    theta = np.linspace(0.2, np.pi - 0.2, n_sensors)
    spos = (
        1.2
        * radius
        * np.c_[np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
    )
    sori = spos / np.linalg.norm(spos, axis=1, keepdims=True)
    squids = os.path.join(out, f"{tag}.squids")
    with open(squids, "w") as fid:
        for p, o in zip(spos, sori):
            fid.write("%g %g %g %g %g %g\n" % (*p, *o))

    ref = np.array(
        [
            [
                sarvas(dipoles[d, :3], dipoles[d, 3:], spos[s], sori[s])
                for d in range(len(dipoles))
            ]
            for s in range(n_sensors)
        ]
    )

    conductivity = (0.3,) if n_layers == 1 else (0.3, 0.006, 0.3)
    meshes = [icosphere(2, radius * f) for f in (1.0, 1.1, 1.2)][:n_layers]

    geom = om.make_nested_geometry(meshes, conductivity=conductivity)
    integrator = om.Integrator(3, 0, 0.005)
    headmat = om.HeadMat(geom, integrator)
    hminv = headmat.inverse()
    dip = om.Matrix(np.asfortranarray(dipoles))
    dsm = om.DipSourceMat(geom, dip, integrator, "Brain")
    sensors = om.Sensors(squids)
    h2mm = om.Head2MEGMat(geom, sensors)
    ds2mm = om.DipSource2MEGMat(dip, sensors)
    gain = om.GainMEG(hminv, dsm, h2mm, ds2mm).array()

    arrays = {
        "ref": ref,
        "dipoles": dipoles,
        "sensor_pos": spos,
        "sensor_ori": sori,
        "HeadMat_flat": np.asarray(headmat.array_flat()),
        "HeadMatInv_flat": np.asarray(hminv.array_flat()),
        "DipSourceMat": dsm.array(),
        "Head2MEGMat": h2mm.array(),
        "DipSource2MEGMat": ds2mm.array(),
        "gain": gain,
    }
    for i, (v, t) in enumerate(meshes):
        arrays[f"mesh{i}_verts"] = v
        arrays[f"mesh{i}_tris"] = t.astype(np.int64)

    # Read the geometry and sensors back out of OpenMEEG rather than trusting
    # what we handed in.  HeadMat and Head2MEGMat have so far been bit-identical
    # on every runner, so if a failing one diverges it must do so either in the
    # assembly or in the inputs as OpenMEEG actually resolved them -- domain
    # ordering, conductivity assignment, vertex ordering, sensor parsing.
    # NB geom.meshes() is an opaque SwigPyObject and Domain.name() hands back a
    # raw std::string*, so go via the accessors that actually work: the known
    # mesh names, the global vertex list, and conductivities in domain order.
    geom_info = {}
    try:
        arrays["geom_vertices"] = np.array(
            [v.array() for v in geom.vertices()], np.float64
        )
        names = ["Cortex"] if n_layers == 1 else ["Cortex", "Skull", "Head"]
        for i, name in enumerate(names):
            mesh = geom.mesh(name)
            gt = np.array(
                [mesh.triangle(t).array() for t in mesh.triangles()], np.int64
            )
            arrays[f"geom_{name}_tris"] = gt
            geom_info[name] = {"n_tris": int(gt.shape[0])}
        geom_info["domain_conductivities"] = [
            float(d.conductivity()) for d in geom.domains()
        ]
        geom_info["n_domains"] = int(len(geom.domains()))
        geom_info["n_vertices"] = int(len(geom.vertices()))
        geom_info["is_nested"] = bool(geom.is_nested())
        geom_info["nb_current_barrier_triangles"] = int(
            geom.nb_current_barrier_triangles()
        )
    except Exception as exc:
        geom_info["error"] = f"{type(exc).__name__}: {exc}"
    try:
        arrays["sensors_pos_readback"] = sensors.getPositions().array()
        arrays["sensors_ori_readback"] = sensors.getOrientations().array()
    except Exception as exc:
        geom_info["sensor_readback_error"] = f"{type(exc).__name__}: {exc}"

    # Split into the physically meaningful pieces.  For radially oriented
    # sensors on a sphere the true secondary field is exactly zero, so RDM is
    # nothing but the spurious BEM secondary term.
    A = unpack_sym(headmat)
    Ainv = unpack_sym(hminv)
    primary = arrays["DipSource2MEGMat"]
    secondary = (arrays["Head2MEGMat"] @ Ainv) @ arrays["DipSourceMat"]
    arrays["primary"] = primary
    arrays["secondary"] = secondary

    rdm, mag = rdm_mag(gain, ref)
    ev = np.abs(np.linalg.eigvalsh(A))
    cond = float(ev.max() / ev.min())
    summary = {
        "n_layers": n_layers,
        "rdm": rdm.tolist(),
        "mag_minus_1": (mag - 1.0).tolist(),
        "max_rdm": float(rdm.max()),
        "passes_0.02": bool(rdm.max() < 0.02),
        "rdm_primary_only": rdm_mag(primary, ref)[0].tolist(),
        "norm_primary": float(np.linalg.norm(primary)),
        "norm_secondary": float(np.linalg.norm(secondary)),
        "secondary_over_primary": float(
            np.linalg.norm(secondary) / np.linalg.norm(primary)
        ),
        "cond_HeadMat": cond,
        "eig_min": float(ev.min()),
        "eig_max": float(ev.max()),
        "inverse_residual": float(
            np.linalg.norm(A @ Ainv - np.eye(A.shape[0])) / A.shape[0]
        ),
        "geometry": geom_info,
        "inverse_thread_sweep": inverse_thread_sweep(headmat, A, cond),
        "digests": {k: digest(v) for k, v in arrays.items()},
    }

    # Re-assemble the operators from scratch a few times in the same process.
    # These are OpenMP-parallel and use no BLAS, so they must be bit-identical;
    # anything else means a race or uninitialized memory in the assembly, which
    # is now one of the few hypotheses still standing.
    # Record the *magnitude* of any spread, not just a boolean: DipSourceMat is
    # already known to vary by 1 ULP (~1.2e-16 relative, a few dozen entries)
    # because its OpenMP accumulation is not associative, and that is benign.
    # Only a spread far above rounding would be a real race worth chasing.
    builders = {
        "HeadMat": lambda: np.asarray(om.HeadMat(geom, integrator).array_flat()),
        "Head2MEGMat": lambda: om.Head2MEGMat(geom, sensors).array(),
        "DipSourceMat": lambda: om.DipSourceMat(geom, dip, integrator, "Brain").array(),
    }
    reassembly = {}
    for name, build in builders.items():
        first = build()
        scale = np.abs(first).max() or 1.0
        worst, n_diff = 0.0, 0
        for _ in range(2):
            other = build()
            d = np.abs(other - first)
            worst = max(worst, float(d.max()))
            n_diff = max(n_diff, int((d > 0).sum()))
        reassembly[name] = {
            "max_abs_dev": worst,
            "max_rel_dev": worst / scale,
            "n_entries_differing": n_diff,
            "size": int(first.size),
            "beyond_rounding": bool(worst / scale > 1e-12),
        }
    summary["reassembly"] = reassembly

    # Repeat in-process: catches run-to-run nondeterminism (threads, ASLR).
    repeats = []
    for _ in range(3):
        g2 = om.GainMEG(
            om.HeadMat(geom, integrator).inverse(),
            om.DipSourceMat(geom, dip, integrator, "Brain"),
            om.Head2MEGMat(geom, sensors),
            om.DipSource2MEGMat(dip, sensors),
        ).array()
        repeats.append(float(rdm_mag(g2, ref)[0].max()))
    summary["max_rdm_repeats"] = repeats

    np.savez_compressed(os.path.join(out, f"{tag}.npz"), **arrays)
    return summary


def main():
    """Dump both models to <output_dir> and print a human-readable summary."""
    out = sys.argv[1] if len(sys.argv) > 1 else "sarvas_dump"
    label = sys.argv[2] if len(sys.argv) > 2 else "dump"
    os.makedirs(out, exist_ok=True)

    report = {"label": label, "ci_failing_rdm_reference": CI_FAILING_RDM}
    try:
        report["environment"] = environment()
    except Exception:
        report["environment_error"] = traceback.format_exc()

    for n_layers in (1, 3):
        tag = f"{label}_{n_layers}layer"
        try:
            report[f"{n_layers}layer"] = compute(n_layers, out, tag)
        except Exception:
            report[f"{n_layers}layer_error"] = traceback.format_exc()

    path = os.path.join(out, f"{label}_summary.json")
    with open(path, "w") as fid:
        json.dump(report, fid, indent=2, default=str)

    print(f"===== sarvas dump [{label}] =====")
    env = report.get("environment", {})
    print(
        f"  python     : {env.get('python_implementation')} "
        f"{sys.version.split()[0]}  numpy {env.get('numpy_version')}"
    )
    print(f"  cpu        : {env.get('cpu_model')}")
    for lib, info in (env.get("openblas") or {}).items():
        # "scipy_openblas" is numpy's bundled copy; anything else is the one
        # OpenMEEG links, which is the one this investigation is about.
        whose = "numpy" if "scipy_openblas" in os.path.basename(lib) else "OpenMEEG"
        print(f"  openblas   : {os.path.basename(lib)}  [{whose}]")
        for k, v in info.items():
            print(f"      {k}: {v}")
    for n_layers in (1, 3):
        s = report.get(f"{n_layers}layer")
        if not s:
            print(f"  {n_layers}-layer: FAILED, see {label}_summary.json")
            continue
        print(
            f"  {n_layers}-layer max RDM = {s['max_rdm']:.6f} "
            f"({'pass' if s['passes_0.02'] else 'FAIL'} vs 0.02), "
            f"repeats={s['max_rdm_repeats']}"
        )
        print(f"      rdm          = {s['rdm']}")
        print(
            f"      cond={s['cond_HeadMat']:.3e} "
            f"inv_resid={s['inverse_residual']:.3e} "
            f"|sec|/|prim|={s['secondary_over_primary']:.6f}"
        )
        print(f"      HeadMat sha  = {s['digests']['HeadMat_flat'][:32]}")
        print(f"      gain    sha  = {s['digests']['gain'][:32]}")
        for name, r in (s.get("reassembly") or {}).items():
            flag = "  *** BEYOND ROUNDING ***" if r["beyond_rounding"] else ""
            print(
                f"      re-assembly {name:<16} max_rel_dev={r['max_rel_dev']:.3e} "
                f"({r['n_entries_differing']}/{r['size']} entries){flag}"
            )
        sweep = s.get("inverse_thread_sweep") or {}
        cells = [
            f"{k}thr:{v['residual']:.1e}"
            for k, v in sweep.items()
            if isinstance(v, dict) and "residual" in v
        ]
        if cells:
            flag = "  *** WRONG INVERSE ***" if sweep.get("any_broken") else ""
            print(f"      inverse vs threads: {'  '.join(cells)}{flag}")
        elif sweep.get("error"):
            print(f"      inverse thread sweep unavailable: {sweep['error']}")
        g = s.get("geometry") or {}
        if g.get("domain_conductivities"):
            print(
                f"      sigma={g['domain_conductivities']} "
                f"nverts={g.get('n_vertices')} nested={g.get('is_nested')} "
                f"barrier_tris={g.get('nb_current_barrier_triangles')}"
            )
        elif g.get("error"):
            print(f"      geometry readback FAILED: {g['error']}")
    print(f"  wrote {path}")


if __name__ == "__main__":
    try:
        main()
    except Exception:  # never let diagnostics fail the job
        traceback.print_exc()
