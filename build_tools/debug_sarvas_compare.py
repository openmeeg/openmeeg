#!/usr/bin/env python
"""Compare Sarvas dumps from several runners and localize the first difference.

Usage:  python debug_sarvas_compare.py <dir> [<dir> ...]

Each <dir> is an unpacked ``sarvas-dump-N`` artifact.  Prints an overview table
(CPU / OpenBLAS kernel / max RDM / pass-fail) and then, for every pair of runs
that disagree, walks the operators in dependency order and reports the *first*
one that differs -- which is the one worth debugging.
"""

from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

# Dependency order: each entry depends only on those above it, so the first
# mismatch in this list is where the divergence is introduced.
PIPELINE = [
    "mesh0_verts",
    "mesh1_verts",
    "mesh2_verts",
    "mesh0_tris",
    "mesh1_tris",
    "mesh2_tris",
    "sensor_pos",
    "sensor_ori",
    "dipoles",
    "ref",
    # what OpenMEEG actually resolved the inputs to, before any assembly
    "geom_vertices",
    "geom_Cortex_tris",
    "geom_Skull_tris",
    "geom_Head_tris",
    "sensors_pos_readback",
    "sensors_ori_readback",
    "HeadMat_flat",
    "HeadMatInv_flat",
    "DipSourceMat",
    "Head2MEGMat",
    "DipSource2MEGMat",
    "primary",
    "secondary",
    "gain",
]


def load(dirs):
    """Load every ``*_summary.json`` found under the given artifact dirs."""
    runs = []
    for d in dirs:
        for path in sorted(glob.glob(os.path.join(d, "*_summary.json"))):
            with open(path) as fid:
                rep = json.load(fid)
            label = rep.get("label", "?")
            runs.append(
                {
                    "dir": d,
                    "name": f"{os.path.basename(os.path.normpath(d))}/{label}",
                    "label": label,
                    "report": rep,
                    "npz": {
                        n: os.path.join(d, f"{label}_{n}layer.npz") for n in (1, 3)
                    },
                }
            )
    return runs


def overview(runs):
    """Print the per-run table of CPU, OpenBLAS kernel and 3-layer result."""
    print(
        f"{'run':<28} {'cpu':<34} {'openblas kernel':<16} {'3L maxRDM':>10} {'ok':>5}"
        f" {'invResid':>9}"
    )
    print("-" * 108)
    for r in runs:
        env = r["report"].get("environment", {})
        cpu = (env.get("cpu_model") or "?").split(":")[-1].strip()[:33]
        # Report the kernel of the OpenBLAS *OpenMEEG* links, not numpy's
        # bundled scipy-openblas: different version, possibly different
        # dispatch on the same CPU, and irrelevant to the BEM assembly.
        kernel = "?"
        for lib, info in (env.get("openblas") or {}).items():
            if not isinstance(info, dict) or not info.get("openblas_get_corename"):
                continue
            core = info["openblas_get_corename"]
            if "scipy_openblas" in os.path.basename(lib):
                if kernel == "?":
                    kernel = f"({core})"  # parenthesised = numpy's, a fallback
            else:
                kernel = core
                break
        s3 = r["report"].get("3layer") or {}
        rdm = s3.get("max_rdm")
        ok = s3.get("passes_0.02")
        # A large inverse residual would mean OpenBLAS's DSPTRF/DSPTRI returned
        # a silently wrong factorization (Info==0 but garbage) -- the arch- and
        # thread-specific wrong-result class of OpenBLAS#3044.  Since RDM
        # responds to HeadMatInv linearly with gain 1.0, an inverse wrong by
        # ~4x is exactly what would take 0.0061 to 0.0245.
        resid = s3.get("inverse_residual")
        print(
            f"{r['name']:<28} {cpu:<34} {kernel:<16} "
            f"{rdm if rdm is None else format(rdm, '10.6f')} "
            f"{'PASS' if ok else 'FAIL':>5}"
            f" {'' if resid is None else format(resid, '9.2e')}"
            f"{'  <-- INVERSE IS WRONG' if (resid or 0) > 1e-10 else ''}"
        )
        # A re-assembly spread far above rounding would mean a race in the
        # OpenMP assembly; DipSourceMat is known to vary by ~1 ULP benignly.
        for name, ra in (s3.get("reassembly") or {}).items():
            if ra.get("beyond_rounding"):
                print(
                    f"    !! {name} re-assembly spread {ra['max_rel_dev']:.3e} "
                    f"({ra['n_entries_differing']}/{ra['size']} entries)"
                )
        sig = (s3.get("geometry") or {}).get("domain_conductivities")
        if sig is not None and sorted(sig) != sorted([0.0, 0.006, 0.3, 0.3]):
            print(f"    !! unexpected domain conductivities: {sig}")
    print()


# OpenBLAS 0.3.20 (what ubuntu-22.04 ships, and what OpenMEEG links in CI) can
# dispatch x86_64 to these, per DYNAMIC_CORE in its Makefile.system.  Only the
# last four are plausible on Azure's server fleet; the rest are pre-2013 parts.
# Tracking kernels rather than CPU model numbers is what matters, because the
# kernel is what would actually differ numerically -- cf.
# https://github.com/OpenMathLib/OpenBLAS/issues/3044
OPENBLAS_0320_X86_KERNELS = [
    "Haswell",
    "Zen",
    "SkylakeX",
    "CooperLake",
]


def coverage(runs):
    """Report which CPUs and which OpenBLAS kernels have been sampled so far."""
    cpus, om_kernels, np_kernels = {}, {}, {}
    for r in runs:
        env = r["report"].get("environment", {})
        cpu = (env.get("cpu_model") or "?").split(":")[-1].strip()
        s3 = r["report"].get("3layer") or {}
        verdict = "pass" if s3.get("passes_0.02") else "FAIL"
        cpus.setdefault(cpu, set()).add(verdict)
        for lib, info in (env.get("openblas") or {}).items():
            if not isinstance(info, dict) or not info.get("openblas_get_corename"):
                continue
            target = (
                np_kernels if "scipy_openblas" in os.path.basename(lib) else om_kernels
            )
            target.setdefault(info["openblas_get_corename"], set()).add(verdict)

    print("=== coverage ===")
    print(f"CPU models seen ({len(cpus)}):")
    for cpu, v in sorted(cpus.items()):
        print(f"  {'FAIL' if 'FAIL' in v else 'pass':<5} {cpu}")
    print("\nOpenBLAS kernels dispatched, as linked by OpenMEEG:")
    if not om_kernels:
        print("  (none recorded -- dump predates probing OpenMEEG's own libopenblas)")
    for k in OPENBLAS_0320_X86_KERNELS:
        if k in om_kernels:
            print(f"  [x] {k:<12} {'FAIL' if 'FAIL' in om_kernels[k] else 'pass'}")
        else:
            print(f"  [ ] {k:<12} not yet sampled")
    for k in sorted(set(om_kernels) - set(OPENBLAS_0320_X86_KERNELS)):
        print(f"  [x] {k:<12} (unexpected for 0.3.20)")
    if np_kernels:
        print(
            f"\n(numpy's bundled scipy-openblas dispatched: {sorted(np_kernels)} "
            f"-- different version, not the one under investigation)"
        )
    print()


def first_difference(a_path, b_path, name_a, name_b):
    """Report the first operator (in dependency order) that differs."""
    if not (os.path.exists(a_path) and os.path.exists(b_path)):
        print(f"  (missing npz for {name_a} or {name_b})")
        return
    A, B = np.load(a_path), np.load(b_path)
    keys = [k for k in PIPELINE if k in A.files and k in B.files]
    keys += [k for k in A.files if k not in PIPELINE and k in B.files]
    found = False
    for k in keys:
        x, y = A[k], B[k]
        if x.shape != y.shape:
            print(f"  {k:<18} SHAPE {x.shape} vs {y.shape}")
            found = True
            break
        if np.array_equal(x, y):
            continue
        denom = np.abs(x).max() or 1.0
        absd = np.abs(x - y).max()
        print(
            f"  first difference at {k!r}: "
            f"max|d| = {absd:.3e}  rel = {absd / denom:.3e}  "
            f"({'ROUNDING' if absd / denom < 1e-12 else 'STRUCTURAL'})"
        )
        found = True
        break
    if not found:
        print("  all arrays bit-identical")


def main():
    """Compare all dumps passed on the command line."""
    dirs = sys.argv[1:]
    if not dirs:
        print(__doc__)
        raise SystemExit(1)
    runs = load(dirs)
    if not runs:
        print("no *_summary.json found")
        raise SystemExit(1)
    overview(runs)
    coverage(runs)

    def maxrdm(r):
        return (r["report"].get("3layer") or {}).get("max_rdm")

    # Group runs by their 3-layer result; anything that disagrees is worth a diff.
    base = runs[0]
    for other in runs[1:]:
        ra, rb = maxrdm(base), maxrdm(other)
        if ra is None or rb is None:
            continue
        same = abs(ra - rb) < 1e-12
        print(
            f"{base['name']}  vs  {other['name']}   "
            f"({ra:.6f} vs {rb:.6f}, {'same' if same else 'DIFFERENT'})"
        )
        if not same:
            first_difference(
                base["npz"][3], other["npz"][3], base["name"], other["name"]
            )
        print()


if __name__ == "__main__":
    main()
