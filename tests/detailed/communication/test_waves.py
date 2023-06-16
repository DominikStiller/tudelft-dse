from unittest import TestCase

from numpy.testing import assert_allclose

from dse.detailed.communication.waves import (
    wavelength,
    calculate_free_space_loss,
    calculate_travel_time,
)


class TestWaves(TestCase):
    def test_wavelength(self):
        assert_allclose(wavelength(300e6), 1, rtol=1e-3)
        assert_allclose(wavelength(3e9), 0.1, rtol=1e-3)

    def test_calculate_free_space_loss(self):
        assert_allclose(calculate_free_space_loss(float("inf"), 42e6), 0)
        # From AE2111-II slides
        assert_allclose(calculate_free_space_loss(2.755e6, 2.5e9), 1.2e-17, rtol=1e-3)

    def test_calculate_travel_time(self):
        assert_allclose(calculate_travel_time(0), 0)
        # One light year
        assert_allclose(calculate_travel_time(9.4607e15), 365.25 * 86400, rtol=1e-5)
