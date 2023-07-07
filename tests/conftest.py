import pytest
import tempfile

@pytest.fixture(scope="session")
def temporary_folder():
    """
    Create a temporary folder that will be used as output directory.
    """
    fn = tempfile.TemporaryDirectory()
    # yiled temporary folder name to be used as Path object
    yield fn.name