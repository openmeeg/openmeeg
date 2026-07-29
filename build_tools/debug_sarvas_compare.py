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
        f"{'run':<28} {'cpu':<34} {'openblas kernel':<16} {'3L maxRDM':>10} {'ok':>4}"
    )
    print("-" * 98)
    for r in runs:
        env = r["report"].get("environment", {})
        cpu = (env.get("cpu_model") or "?").split(":")[-1].strip()[:33]
        kernel = "?"
        for info in (env.get("openblas") or {}).values():
            if isinstance(info, dict) and info.get("openblas_get_corename"):
                kernel = info["openblas_get_corename"]
                break
        s3 = r["report"].get("3layer") or {}
        rdm = s3.get("max_rdm")
        ok = s3.get("passes_0.02")
        print(
            f"{r['name']:<28} {cpu:<34} {kernel:<16} "
            f"{rdm if rdm is None else format(rdm, '10.6f')} "
            f"{'PASS' if ok else 'FAIL':>4}"
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
