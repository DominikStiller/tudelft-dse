from unittest import TestCase

from numpy.testing import assert_allclose

from dse.detailed.communication import from_db, to_db


class TestCommunication(TestCase):
    def test_db(self):
        for x in [1, 14, 0.53]:
            assert_allclose(x, from_db(to_db(x)))

        for x in [0, 14, -3.5]:
            assert_allclose(x, to_db(from_db(x)))
