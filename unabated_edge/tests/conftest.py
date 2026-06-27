import pytest
from unabated_edge import storage

@pytest.fixture(autouse=True)
def _clean_storage_buffer():
    storage._BUFFER.clear()
    yield
    storage._BUFFER.clear()
