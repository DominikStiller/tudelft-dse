from unittest import TestCase

from numpy.testing import assert_allclose


class TestExample(TestCase):
    def test_function(self):
        self.assertEquals(1, 1)
        assert_allclose(0, 0.01)
