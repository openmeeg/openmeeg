"""Regression test for https://github.com/openmeeg/openmeeg/issues/818.

Inverting ``HeadMat`` for a large head model (there, dimension 59420) used to
segfault: with a 32-bit-integer BLAS/LAPACK the packed symmetric factorization
(``DSPTRF``/``DSPTRI``) overflows its internal indices once the matrix dimension
exceeds ~46340, and reads out of bounds. The fix is to build/link against an
ILP64 (64-bit integer) OpenBLAS -- i.e. ``scipy-openblas64``.

Actually reproducing the crash requires a matrix past that boundary: the packed
storage alone is >8.6 GB and the O(n^3) inversion takes hours, so it is not a
practical unit test. Instead we assert the *fix* is in place: that OpenMEEG's
own (vendored) OpenBLAS is the 64-bit-integer build. If a build regresses to the
LP64 ``scipy-openblas32`` this fails; on builds that link an external/system
OpenBLAS (no vendored library, e.g. the source builds) it skips.
"""

import os

import pytest

import openmeeg  # noqa: F401  -- importing loads the extension and its BLAS


def test_issue_818_openmeeg_uses_ilp64_openblas():
    threadpoolctl = pytest.importorskip("threadpoolctl")

    def _norm(path):
        return path.replace("\\", "/").lower()

    infos = threadpoolctl.threadpool_info()
    # Restrict to OpenMEEG's *own* BLAS, which is vendored next to the
    # extension module (e.g. site-packages/openmeeg.libs/, openmeeg/.dylibs/).
    # This avoids mistaking, say, NumPy's OpenBLAS for OpenMEEG's.
    om_blas = [
        info
        for info in infos
        if info.get("internal_api") == "openblas"
        and "openmeeg" in _norm(info.get("filepath", ""))
    ]
    if not om_blas:
        pytest.skip(
            "OpenMEEG is not linked against a vendored OpenBLAS (e.g. a source "
            "build against a system OpenBLAS); the ILP64 check does not apply."
        )

    for info in om_blas:
        name = os.path.basename(_norm(info["filepath"]))
        # scipy-openblas64 is named libscipy_openblas64_ and its symbols carry a
        # 64_ suffix; scipy-openblas32 (LP64) has neither.
        assert "openblas64" in name or "64_" in name, (
            f"OpenMEEG links a non-ILP64 OpenBLAS ({name!r}); issue #818 "
            "requires scipy-openblas64 (64-bit BLAS integers)."
        )
