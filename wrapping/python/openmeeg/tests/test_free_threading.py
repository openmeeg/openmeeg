import concurrent.futures as cf
import os.path as op
import sys
import sysconfig

import numpy as np
import pytest
from numpy.testing import assert_allclose

import openmeeg as om
from openmeeg._openmeeg_wrapper import Vector

free_threaded = bool(sysconfig.get_config_var("Py_GIL_DISABLED"))


@pytest.mark.skipif(not free_threaded, reason="requires a free-threaded build")
def test_gil_not_reenabled():
    """The extension must not make CPython fall back to enabling the GIL."""
    # Importing a module without Py_MOD_GIL_NOT_USED re-enables the GIL process-wide
    # rather than failing, so this is the only thing that catches a missing -nogil.
    assert not sys._is_gil_enabled()


@pytest.mark.skipif(not free_threaded, reason="requires a free-threaded build")
def test_concurrent_geometry_load(data_path):
    """Load geometries concurrently: gh-866 made the i/o registries safe for this."""
    geom_file = op.join(data_path, "Head1", "Head1.geom")
    cond_file = op.join(data_path, "Head1", "Head1.cond")

    def load():
        return om.read_geometry(geom_file, cond_file).nb_parameters()

    with cf.ThreadPoolExecutor(max_workers=8) as pool:
        results = [f.result() for f in [pool.submit(load) for _ in range(32)]]

    assert len(set(results)) == 1, results
    assert not sys._is_gil_enabled()


@pytest.mark.skipif(not free_threaded, reason="requires a free-threaded build")
def test_concurrent_typemap_conversion():
    """Exercise the array typemaps from several threads at once."""
    arrays = [np.full(64, float(i), order="F") for i in range(32)]

    def roundtrip(arr):
        # Vector.array() returns a view, so copy before the Vector goes away
        vec = Vector(arr)
        return vec.array().copy()

    with cf.ThreadPoolExecutor(max_workers=8) as pool:
        out = list(pool.map(roundtrip, arrays))

    for got, want in zip(out, arrays):
        assert_allclose(got, want)
