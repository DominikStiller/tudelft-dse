from typing import Self

import pandas as pd
import scipy.optimize
from matplotlib import pyplot as plt
from tabulate import tabulate, SEPARATING_LINE

from dse.detailed.communication import from_db, to_db, k
from dse.detailed.communication.waves import (
    calculate_free_space_loss,
    wavelength,
    calculate_travel_time,
)
from dse.plotting import format_plot


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


def calculate_system_noise_temp(line_loss: float, noise_figure: float, temperature_antenna: float):
    temperature_line = calculate_line_temp(line_loss)
    temperature_receiver = calculate_receiver_temp(line_loss, noise_figure)
    return temperature_antenna + temperature_line + temperature_receiver


def calculate_uplink_data_rate():
    data_rate_voice = 8e3  # CCSDS 766.2-B-1, 4.5.1, use G.729 with 8 kb/s
    data_rate_biomed = 0.2e3  # 50 * 4 * 1
    data_rate_telemetry = 8e3  # 200 * 4 * 10
    factor = 2
    return factor * (data_rate_voice + data_rate_biomed + data_rate_telemetry)


def calculate_downlink_data_rate():
    data_rate_voice = 8e3  # CCSDS 766.2-B-1, 4.5.1, use G.729 with 8 kb/s
    data_rate_mission = 10e3
    factor = 2
    return factor * (data_rate_voice + data_rate_mission)


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
    elif coding == "ldpc-1/2":
        # From CCSDS 130.1-G-3, Figure 3-5
        return 1
    else:
        raise "Unsupported coding scheme"


def coding_datarate_factor(coding: str):
    if coding == "uncoded":
        return 1
    elif "1/2" in coding:
        return 1 / 2
    else:
        raise "Unsupported coding scheme"


class Link:
    def __init__(
        self,
        name: str,
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
        data_rate: float,
        bit_error_rate: float,
        coding: str,
        line_loss: float = None,
        noise_figure: float = None,
        temperature_antenna: float = None,
        temperature_system_noise: float = None,
        margin_target: float = 3,
        latency_factor: float = 1,
    ):
        """
        Initialize class.

        Args:
            frequency:
            tx_power:
            tx_loss: [dB]
            tx_gain: [dBi]
            tx_pointing_loss: negative [dBi]
            tx_rx_distance:
            loss_environment: negative [dBi]
            rx_pointing_loss: negative [dBi]
            rx_gain: [dBi]
            rx_loss: [dB]
            data_rate:
            bit_error_rate:
            coding:
            line_loss: [-]
            noise_figure: [-]
            temperature_antenna:
            temperature_system_noise:
            margin_target: [dB]
            latency_factor:
        """
        self.name = name
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
        self.data_rate = data_rate
        self.bit_error_rate = bit_error_rate
        self.coding = coding
        self.margin_target = margin_target
        self.latency_factor = latency_factor

        if temperature_system_noise is not None:
            assert line_loss is None
            assert noise_figure is None
            assert temperature_antenna is None
            self.temperature_system_noise = temperature_system_noise
        else:
            self.line_loss = line_loss
            self.noise_figure = noise_figure
            self.temperature_antenna = temperature_antenna
            self.temperature_system_noise = calculate_system_noise_temp(
                line_loss, noise_figure, temperature_antenna
            )

        self.budget = None
        self.snr = None
        self.max_data_rate = None

        self._calculate_max_data_rate()
        self._calculate_budget(self.data_rate)

    def _calculate_budget(self, data_rate):
        loss_free_space = to_db(calculate_free_space_loss(self.tx_rx_distance, self.frequency))
        required_datarate = data_rate / coding_datarate_factor(self.coding)

        self.budget = pd.DataFrame(
            [
                ("Tx power", to_db(self.tx_power), "dBW"),
                ("Tx line loss", self.tx_loss, "dB"),
                ("Tx antenna gain", self.tx_gain, "dBi"),
                ("Tx antenna pointing loss", self.tx_pointing_loss, "dB"),
                ("Free space loss", loss_free_space, "dB"),
                ("Environment loss", self.loss_environment, "dB"),
                ("Rx antenna gain", self.rx_gain, "dBi"),
                ("Rx antenna pointing loss", self.rx_pointing_loss, "dB"),
                ("Rx line loss", self.rx_loss, "dB"),
                ("Required data rate", to_db(1 / required_datarate), "dB(bit/s)"),
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

        return snr_margin

    def _calculate_max_data_rate(self):
        def fn(data_rate):
            return self._calculate_budget(data_rate) - self.margin_target

        res = scipy.optimize.least_squares(fn, self.data_rate)
        self.max_data_rate = res.x[0]

    def print_configuration(self):
        print(self.name.upper())
        print(f"Transmitter power: {self.tx_power:.0f} W")
        print(f"Frequency: {self.frequency / 1e6:.1f} MHz")
        print(f"System noise temperature: {self.temperature_system_noise:.0f} K")
        print(f"Wavelength: {wavelength(self.frequency):.3g} m")
        print(f"1/4 Wavelength: {wavelength(self.frequency)/4:.3g} m")

        g_over_t = self.rx_gain - to_db(self.temperature_system_noise)
        print(f"G/T: {g_over_t:.1f} dB/K")

        print(f"Nominal data rate: {self.data_rate/1e3:.1f} kb/s")
        print(
            f"Maximum data rate (with {self.margin_target:.0f} dB margin): {self.max_data_rate/1e3:.1f} kb/s"
        )

        latency = self.latency_factor * calculate_travel_time(self.tx_rx_distance)
        print(f"Latency: {latency * 1e3:.0f} ms")

    def print_table(self):
        budget = list(self.budget.itertuples(index=False, name=None))
        snr = list(self.snr.itertuples(index=False, name=None))
        print(tabulate(budget + [SEPARATING_LINE] + snr, numalign="right", floatfmt="+.3f"))

    @classmethod
    def _make_plotting_df(cls):
        return pd.DataFrame(
            [
                ("Tx power", ["Tx power"], 0),
                ("Antenna", ["Tx line loss"], 0.15),
                ("EIRP", ["Tx antenna gain", "Tx antenna pointing loss"], 0.3),
                ("Received radiation", ["Free space loss", "Environment loss"], 0.7),
                ("Rx power", ["Rx antenna gain", "Rx antenna pointing loss"], 0.85),
                ("Rx power 2", ["Rx line loss"], 1),
            ],
            columns=["loc_name", "components", "pos"],
        )

    def plot(self, ax: plt.Axes = None):
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 4))

        df = self._make_plotting_df()

        def get_component_value(df: pd.DataFrame, component: str):
            return df[df["Component"] == component]["Value [dB]"].iloc[0]

        def get_gain(row):
            gain = 0
            for component in row["components"]:
                gain += get_component_value(self.budget, component)
            return gain

        df["gain"] = df.apply(get_gain, axis=1)
        df["cum_gain"] = df["gain"].cumsum()

        power_required = (
            get_component_value(self.snr, "Required Eb/N0")
            + self.margin_target
            - get_component_value(self.budget, "Required data rate")
            - get_component_value(self.budget, "Boltzmann constant")
            - get_component_value(self.budget, "System noise temperature")
        )

        for pos in df["pos"]:
            ax.axvline(pos, color="black", alpha=0.5, ls=":")

        ax.plot(df["pos"], df["cum_gain"], c="black")
        ax.hlines(
            power_required,
            0.95,
            1.05,
            ls="--",
            color="black",
            label=f"Required power (with {self.margin_target} dB margin)",
        )

    @classmethod
    def plot_multiple(cls, links: list[Self]):
        _, axs = plt.subplots(len(links), figsize=(10, 2 * len(links)), sharex="all", sharey="all")

        df = cls._make_plotting_df()

        for ax, link in zip(axs, links):
            link.plot(ax)
            ax.set_title(link.name)
            ax.set_ylabel("Power [dB]")

        axs[0].legend(loc="upper right")
        axs[-1].set_xticks(df["pos"], df["loc_name"])

        format_plot()
        plt.show()


if __name__ == "__main__":
    relay_gain_uhf = 0
    aircraft_gain_uhf = 4
    aircraft_gain_hf = 0
    base_gain_hf = 0

    tx_loss = -1
    rx_loss = -2

    relay_tx_rx_distance = 8000e3
    skywave_tx_rx_distance = 933e3

    budget_uplink_relay = Link(
        name="Uplink (relay)",
        frequency=440e6,
        tx_power=55,
        tx_loss=tx_loss,
        tx_gain=aircraft_gain_uhf,
        tx_pointing_loss=-2,
        tx_rx_distance=relay_tx_rx_distance,
        loss_environment=-0.5,
        rx_pointing_loss=-2,
        rx_gain=relay_gain_uhf,
        rx_loss=rx_loss,
        data_rate=calculate_uplink_data_rate(),
        bit_error_rate=1e-6,
        coding="ldpc-1/2",
        temperature_system_noise=500,
        latency_factor=2,
    )
    budget_uplink_relay.print_configuration()
    budget_uplink_relay.print_table()

    print()
    budget_downlink_relay = Link(
        name="Downlink (relay)",
        frequency=420e6,
        tx_power=100,
        tx_loss=tx_loss,
        tx_gain=relay_gain_uhf,
        tx_pointing_loss=-2,
        tx_rx_distance=relay_tx_rx_distance,
        loss_environment=-0.5,
        rx_pointing_loss=-2,
        rx_gain=aircraft_gain_uhf,
        rx_loss=rx_loss,
        line_loss=from_db(-1),
        noise_figure=from_db(4.9),
        temperature_antenna=54,
        data_rate=calculate_downlink_data_rate(),
        bit_error_rate=1e-6,
        coding="ldpc-1/2",
        latency_factor=2,
    )
    budget_downlink_relay.print_configuration()
    budget_downlink_relay.print_table()

    print()
    budget_uplink_skywave = Link(
        name="Uplink (skywave)",
        frequency=4e6,
        tx_power=1,
        tx_loss=tx_loss,
        tx_gain=aircraft_gain_hf,
        tx_pointing_loss=-2,
        tx_rx_distance=skywave_tx_rx_distance,
        loss_environment=2 * -0.5,
        rx_pointing_loss=-2,
        rx_gain=base_gain_hf,
        rx_loss=rx_loss,
        line_loss=from_db(-1),
        noise_figure=from_db(4.9),
        temperature_antenna=155,
        data_rate=calculate_uplink_data_rate(),
        bit_error_rate=1e-6,
        coding="ldpc-1/2",
    )
    budget_uplink_skywave.print_configuration()
    budget_uplink_skywave.print_table()

    print()
    budget_downlink_skywave = Link(
        name="Downlink (skywave)",
        frequency=4e6,
        tx_power=1,
        tx_loss=tx_loss,
        tx_gain=base_gain_hf,
        tx_pointing_loss=-2,
        tx_rx_distance=skywave_tx_rx_distance,
        loss_environment=2 * -0.5,
        rx_pointing_loss=-2,
        rx_gain=aircraft_gain_hf,
        rx_loss=rx_loss,
        line_loss=from_db(-1),
        noise_figure=from_db(4.9),
        temperature_antenna=54,
        data_rate=calculate_downlink_data_rate(),
        bit_error_rate=1e-6,
        coding="ldpc-1/2",
    )
    budget_downlink_skywave.print_configuration()
    budget_downlink_skywave.print_table()

    Link.plot_multiple([budget_uplink_relay, budget_uplink_skywave])
    # Link.plot_multiple([budget_downlink_relay, budget_downlink_skywave])
