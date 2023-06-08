from unittest import TestCase

from numpy.testing import assert_allclose

from dse.detailed.communication import from_db
from dse.detailed.communication.linkbudget import (
    calculate_line_temp,
    calculate_receiver_temp,
    calculate_system_noise_temp,
    calculate_uplink_data_rate,
    calculate_downlink_data_rate,
    required_snr,
    coding_datarate_factor,
    Link,
)


class TestLinkbudget(TestCase):
    def test_calculate_line_temp(self):
        assert_allclose(calculate_line_temp(1), 0)
        # From Wertz p. 558
        assert_allclose(calculate_line_temp(0.89), 36, rtol=1e-2)

    def test_calculate_receiver_temp(self):
        assert_allclose(calculate_receiver_temp(0.5, 1), 0)
        # From Wertz p. 558
        assert_allclose(calculate_receiver_temp(0.89, 1.1), 33, rtol=1e-1)

    def test_calculate_system_noise_temp(self):
        assert_allclose(calculate_system_noise_temp(1, 1, 0), 0)
        # From Wertz p. 558
        assert_allclose(
            calculate_system_noise_temp(from_db(-0.5), from_db(0.5), 150), 221, rtol=1e-1
        )
        assert_allclose(
            calculate_system_noise_temp(from_db(-0.5), from_db(3.0), 290), 614, rtol=1e-1
        )

    def test_calculate_uplink_data_rate(self):
        assert_allclose(calculate_uplink_data_rate(), 32.4e3)

    def test_calculate_downlink_data_rate(self):
        assert_allclose(calculate_downlink_data_rate(), 36e3)

    def test_required_snr(self):
        assert_allclose(required_snr(1e-6, "bpsk", "uncoded"), 10.5)
        assert_allclose(required_snr(1e-6, "bpsk", "ldpc-1/2"), 1.0)
        assert_allclose(required_snr(1e-6, "8fsk", "uncoded"), 10.0)

    def test_coding_datarate_factor(self):
        assert_allclose(coding_datarate_factor("uncoded"), 1)
        assert_allclose(coding_datarate_factor("ldpc-1/2"), 0.5)

    def test_link(self):
        # From AE2111-II slides
        link = Link(
            "Test",
            frequency=2.5e9,
            tx_power=2,
            tx_loss=-0.97,
            tx_gain=19.7,
            tx_pointing_loss=-0.003,
            tx_rx_distance=2.755e6,
            loss_environment=-0.5,
            rx_pointing_loss=-0.12,
            rx_gain=45.8,
            rx_loss=-0.97,
            data_rate=6e6,
            bit_error_rate=1e-6,
            modulation="8fsk",
            coding="uncoded",
            temperature_system_noise=135,
        )

        assert_allclose(link.snr.loc["Received Eb/N0"][0], 36.2, 1e-2)
        assert_allclose(link.snr.loc["Margin"][0], 26.2, 1e-2)
