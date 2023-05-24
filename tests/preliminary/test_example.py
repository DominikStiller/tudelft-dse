from unittest import TestCase

from numpy.testing import assert_allclose


class TestExample(TestCase):
    def test_function(self):
        self.assertEqual(1, 1)
        assert_allclose(1.0000001e9, 1.0000002e9)
