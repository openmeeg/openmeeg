import os
import pytest


@pytest.fixture(scope='session')
def data_path():
    data_path = os.getenv('OPENMEEG_DATA_PATH')
    assert data_path is not None, 'OPENMEEG_DATA_PATH must be set'
    assert os.path.isdir(data_path), 'OPENMEEG_DATA_PATH must exist'
    return data_path
