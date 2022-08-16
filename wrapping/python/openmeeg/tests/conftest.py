import os
import pytest


@pytest.fixture(scope="session")
def data_path():
    data_path = os.getenv("OPENMEEG_DATA_PATH")
    assert data_path is not None, "OPENMEEG_DATA_PATH must be set, got None"
    # deal with MSVC not handling mixed paths like
    # D:/a/openmeeg/openmeeg/data\Head1\Head1.geom
    # but cmake uses a mixed path for the --path arg
    data_path = data_path.replace("/", os.path.sep)
    assert os.path.isdir(data_path), f"OPENMEEG_DATA_PATH does not exist: ${data_path}"
    return data_path
