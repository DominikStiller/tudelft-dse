from unittest import TestCase

from dse.detailed.communication.antenna import from_db, to_db


class TestAntenna(TestCase):
    def test_db(self):
        for x in [1, 14, 0.53]:
            self.assertAlmostEquals(x, from_db(to_db(x)))

        for x in [0, 14, -3.5]:
            self.assertAlmostEquals(x, to_db(from_db(x)))
