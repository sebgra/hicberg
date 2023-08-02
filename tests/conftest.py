import pytest
import tempfile
import numpy
import random as rand

@pytest.fixture(scope="session")
def temporary_folder():
    """
    Create a temporary folder that will be used as output directory.
    """
    fn = tempfile.TemporaryDirectory()
    # yiled temporary folder name to be used as Path object
    yield fn.name


@pytest.fixture
def random():
    rand.seed(0)
    yield numpy.random.seed(0)