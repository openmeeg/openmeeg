import os
from pathlib import Path

import numpy as np
import pytest


def pytest_configure(config):
    """Configure pytest options."""
    config.addinivalue_line("usefixtures", "run_some_parallel")
    config.addinivalue_line("markers", "slow: marks tests as slow")


@pytest.fixture(scope="session")
def data_path():
    data_path = os.getenv("OPENMEEG_DATA_PATH")
    assert data_path is not None, "OPENMEEG_DATA_PATH must be set, got None"
    # deal with MSVC not handling mixed paths like
    # D:/a/openmeeg/openmeeg/data\Head1\Head1.geom
    # but cmake uses a mixed path for the --path arg
    data_path = data_path.replace("/", os.path.sep)
    assert os.path.isdir(data_path), f"OPENMEEG_DATA_PATH does not exist: ${data_path}"
    return Path(data_path)


@pytest.fixture(scope="session")
def run_some_parallel():
    """Run some stuff in parallel."""
    # This is to try to get some NumPy parallelism in the tests.
    a = np.ones((1000, 1000))
    _ = a @ a
