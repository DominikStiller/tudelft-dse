import pandas as pd
from tabulate import tabulate, SEPARATING_LINE

from dse.detailed.communication import from_db, to_db, k
from dse.detailed.communication.waves import calculate_free_space_loss, wavelength


def calculate_line_temp(line_loss: float, T0: float = 290):
    # line_loss [-]
    # T0 [K]
    assert 0 < line_loss <= 1
    return T0 * (1 - line_loss) / line_loss


def calculate_receiver_temp(line_loss: float, noise_figure: float, T0: float = 290):
    # line_loss [-]
    # noise_figure [-]
    # T0 [K]
    assert 0 < line_loss <= 1
    assert noise_figure >= 1
    return T0 * (noise_figure - 1) / line_loss


def calculate_uplink_data_rate():
    data_rate_voice = 128e3  # CCSDS 766.2-B-1, 4.2.1.5
    data_rate_biomed = 33e3
    data_rate_telemetry = 0.8e3
    factor = 1.5
    return factor * (data_rate_voice + data_rate_biomed + data_rate_telemetry)


def required_snr(bit_error_rate: float, coding: str):
    if bit_error_rate != 1e-6:
        raise "Unsupported BER"

    if coding == "uncoded":
        # From CCSDS 130.1-G-3, Figure 3-5
        return 10.5
    elif coding == "turbo-8920-1/2":
        # From CCSDS 130.1-G-3, Figure 7-9
        return 1.1
    elif coding == "conv-7-1/2":
        # From CCSDS 130.1-G-3, Figure 3-5
        return 4.8
    else:
        raise "Unsupported coding scheme"


class LinkBudget:
    def __init__(
        self,
        frequency: float,
        tx_power: float,
        tx_loss: float,
        tx_gain: float,
        tx_pointing_loss: float,
        tx_rx_distance: float,
        loss_environment: float,
        rx_pointing_loss: float,
        rx_gain: float,
        rx_loss: float,
        line_loss: float,
        noise_figure: float,
        temperature_antenna: float,
        data_rate: float,
        bit_error_rate: float,
        coding: str,
    ):
        """
        Initialize class.

        Args:
            frequency:
            tx_power:
            tx_loss:
            tx_gain: [dBi]
            tx_pointing_loss: negative [dBi]
            tx_rx_distance:
            loss_environment: negative [dBi]
            rx_pointing_loss: negative [dBi]
            rx_gain: [dBi]
            rx_loss: [-]
            line_loss: [-]
            noise_figure: [-]
            temperature_antenna:
            data_rate:
            bit_error_rate:
            coding:
        """
        self.frequency = frequency
        self.tx_power = tx_power
        self.tx_loss = tx_loss
        self.tx_gain = tx_gain
        self.tx_pointing_loss = tx_pointing_loss
        self.tx_rx_distance = tx_rx_distance
        self.loss_environment = loss_environment
        self.rx_pointing_loss = rx_pointing_loss
        self.rx_gain = rx_gain
        self.rx_loss = rx_loss
        self.line_loss = line_loss
        self.noise_figure = noise_figure
        self.temperature_antenna = temperature_antenna
        self.data_rate = data_rate
        self.bit_error_rate = bit_error_rate
        self.coding = coding

        self.budget = None
        self.temperature_system_noise = None
        self.snr = None

        self._calculate_budget()

    def _calculate_budget(self):
        loss_free_space = to_db(calculate_free_space_loss(self.tx_rx_distance, self.frequency))

        temperature_line = calculate_line_temp(self.line_loss)
        temperature_receiver = calculate_receiver_temp(self.line_loss, self.noise_figure)
        self.temperature_system_noise = (
            self.temperature_antenna + temperature_line + temperature_receiver
        )
        # self.temperature_system_noise = 500

        self.budget = [
            to_db(self.tx_power),
            to_db(self.tx_loss),
            self.tx_gain,
            self.tx_pointing_loss,
            loss_free_space,
            self.loss_environment,
            self.rx_pointing_loss,
            self.rx_gain,
            self.rx_loss,
        ]
        self.budget = pd.DataFrame(
            [
                ("Tx power", to_db(self.tx_power), "dBW"),
                ("Tx loss factor", to_db(self.tx_loss), "dB"),
                ("Tx antenna pointing loss", self.tx_pointing_loss, "dB"),
                ("Tx antenna gain", self.tx_gain, "dBi"),
                ("Free space loss", loss_free_space, "dB"),
                ("Environment loss", self.loss_environment, "dB"),
                ("Rx antenna gain", self.rx_gain, "dBi"),
                ("Rx antenna pointing loss", self.rx_pointing_loss, "dB"),
                ("Rx loss factor", to_db(self.rx_loss), "dB"),
                ("Required data rate", to_db(1 / self.data_rate), "dB(bit/s)"),
                ("Boltzmann constant", to_db(1 / k), "dB(J/K)"),
                ("System noise temperature", to_db(1 / self.temperature_system_noise), "dB(K)"),
            ],
            columns=["Component", "Value [dB]", "Unit"],
        )

        snr_required = required_snr(self.bit_error_rate, self.coding)
        snr_received = self.budget["Value [dB]"].sum()
        snr_margin = snr_received - snr_required
        self.snr = pd.DataFrame(
            [
                ("Received Eb/N0", snr_received, "dB"),
                ("Required Eb/N0", snr_required, "dB"),
                ("Margin", snr_margin, "dB"),
            ],
            columns=["Component", "Value [dB]", "Unit"],
        )

    def print_configuration(self):
        print(f"Transmitter power: {self.tx_power:.0f} W")
        print(f"Frequency: {self.frequency / 1e6:.1f} MHz")
        print(f"Wavelength: {wavelength(self.frequency):.3g} m")
        print(f"1/4 Wavelength: {wavelength(self.frequency)/4:.3g} m")

        g_over_t = self.rx_gain - to_db(self.temperature_system_noise)
        print(f"G/T: {g_over_t:.1f} dB/K")

    def print_table(self):
        budget = list(self.budget.itertuples(index=False, name=None))
        snr = list(self.snr.itertuples(index=False, name=None))
        print(tabulate(budget + [SEPARATING_LINE] + snr, numalign="right", floatfmt="+.3f"))


if __name__ == "__main__":
    aircraft_gain_uhf = 6
    aircraft_gain_hf = 0

    print("UPLINK (RELAY)")
    budget_uplink_relay = LinkBudget(
        frequency=440e6,
        tx_power=100,
        tx_loss=0.7,
        tx_gain=aircraft_gain_uhf,
        tx_pointing_loss=-2,
        tx_rx_distance=8000e3,
        loss_environment=-0.5,
        rx_pointing_loss=-2,
        rx_gain=0,
        rx_loss=0.7,
        line_loss=from_db(-1),
        noise_figure=from_db(4.9),
        temperature_antenna=155,
        data_rate=calculate_uplink_data_rate(),
        bit_error_rate=1e-6,
        coding="conv-7-1/2",
    )
    budget_uplink_relay.print_configuration()
    budget_uplink_relay.print_table()

    print()
    print("UPLINK (SKYWAVE)")
    budget_uplink_skywave = LinkBudget(
        frequency=10e6,
        tx_power=20,
        tx_loss=0.7,
        tx_gain=aircraft_gain_hf,
        tx_pointing_loss=-2,
        tx_rx_distance=933e3,
        loss_environment=-0.5,
        rx_pointing_loss=-2,
        rx_gain=0,
        rx_loss=0.7,
        line_loss=from_db(-1),
        noise_figure=from_db(4.9),
        temperature_antenna=155,
        data_rate=calculate_uplink_data_rate(),
        bit_error_rate=1e-6,
        coding="conv-7-1/2",
    )
    # TODO check if losses are correct
    budget_uplink_skywave.print_configuration()
    budget_uplink_skywave.print_table()

    # TODO downlink
    # TODO calculate max data rate
