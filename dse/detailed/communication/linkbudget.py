from typing import Self

import pandas as pd
from matplotlib import pyplot as plt
from tabulate import tabulate, SEPARATING_LINE
from termcolor import colored

from dse.detailed.communication import from_db, to_db, k, db_to_dbm, dbm_to_db
from dse.detailed.communication.waves import (
    calculate_free_space_loss,
    wavelength,
    calculate_travel_time,
)
from dse.plotting import format_plot, save_plot


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


def calculate_uplink_data_rate(factor=2):
    data_rate_voice = 8e3  # CCSDS 766.2-B-1, 4.5.1, use G.729 with 8 kb/s
    data_rate_biomed = 0.2e3  # 50 * 4 * 1
    data_rate_telemetry = 8e3  # 200 * 4 * 10
    return factor * (data_rate_voice + data_rate_biomed + data_rate_telemetry)


def calculate_downlink_data_rate(factor=2):
    data_rate_voice = 8e3  # CCSDS 766.2-B-1, 4.5.1, use G.729 with 8 kb/s
    data_rate_mission = 10e3
    return factor * (data_rate_voice + data_rate_mission)


def required_snr(bit_error_rate: float, modulation: str, coding: str):
    if bit_error_rate != 1e-6:
        raise "Unsupported BER"

    if modulation == "bpsk":
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
        elif coding == "ldpc-1/2-gladden":
            # From CCSDS 130.1-G-3, Figure 3-5
            return 2.5
        elif coding == "ldpc-7/8":
            # From CCSDS 130.1-G-3, Figure 8-8
            return 3.8
        elif coding == "turbo-1/6":
            # Only used for MarCO validation
            return 3.5
    elif modulation == "8fsk":
        if coding == "uncoded":
            # From Wertz Figure 13-9
            return 10
    elif modulation == "8bsk":
        if coding == "uncoded":
            # From Sumic 2010
            return 11.5

    raise "Unsupported combination of modulation and coding schemes"


def coding_datarate_factor(coding: str):
    if coding == "uncoded":
        return 1
    elif "1/2" in coding:
        return 1 / 2
    elif "1/6" in coding:
        return 1 / 6
    elif "7/8" in coding:
        return 7 / 8
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
        modulation: str,
        coding: str,
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
            modulation:
            coding:
            noise_figure: [dB]
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
        self.modulation = modulation
        self.coding = coding
        self.margin_target = margin_target
        self.latency_factor = latency_factor

        if temperature_system_noise is not None:
            assert noise_figure is None
            assert temperature_antenna is None
            self.temperature_antenna = None
            self.temperature_system_noise = temperature_system_noise
        else:
            self.temperature_antenna = temperature_antenna
            self.temperature_line = calculate_line_temp(from_db(rx_loss))
            self.temperature_receiver = calculate_receiver_temp(
                from_db(rx_loss), from_db(noise_figure)
            )
            self.temperature_system_noise = calculate_system_noise_temp(
                from_db(rx_loss), from_db(noise_figure), temperature_antenna
            )

        self.budget = None
        self.snr = None
        self.max_data_rate = None
        self.eirp = None
        self.rx_power = None

        self._calculate_budget()
        self._validate_budget()

    def _calculate_budget(self):
        loss_free_space = to_db(calculate_free_space_loss(self.tx_rx_distance, self.frequency))
        required_datarate = self.data_rate  # / coding_datarate_factor(self.coding)

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
                ("Data rate", to_db(1 / required_datarate), "dB(bit/s)"),
                (
                    "Noise spectral density",
                    to_db(1 / (k * self.temperature_system_noise)),
                    "dB(W/Hz)",
                ),
            ],
            columns=["Component", "Value [dB]", "Unit"],
        ).set_index("Component", drop=True)

        self.eirp = to_db(self.tx_power) + self.tx_loss + self.tx_gain
        self.rx_power = self.budget.iloc[:-2]["Value [dB]"].sum()

        snr_required = required_snr(self.bit_error_rate, self.modulation, self.coding)
        snr_received = self.budget["Value [dB]"].sum()
        snr_margin = snr_received - snr_required
        self.snr = pd.DataFrame(
            [
                ("Received Eb/N0", snr_received, "dB"),
                ("Required Eb/N0", snr_required, "dB"),
                ("Margin", snr_margin, "dB"),
            ],
            columns=["Component", "Value [dB]", "Unit"],
        ).set_index("Component", drop=True)

        return snr_margin

    def _validate_budget(self):
        assert self._get_component_value(self.budget, "Tx line loss") <= 0
        assert self._get_component_value(self.budget, "Tx antenna gain") >= 0
        assert self._get_component_value(self.budget, "Tx antenna pointing loss") <= 0
        assert self._get_component_value(self.budget, "Free space loss") <= 0
        assert self._get_component_value(self.budget, "Environment loss") <= 0
        assert self._get_component_value(self.budget, "Rx antenna gain") >= 0
        assert self._get_component_value(self.budget, "Rx antenna pointing loss") <= 0
        assert self._get_component_value(self.budget, "Rx line loss") <= 0

    def print_configuration(self):
        print(self.name.upper())
        print(f"Transmitter power: {self.tx_power:.0f} W")
        print(f"EIRP: {self.eirp:.0f} dBW = {db_to_dbm(self.eirp):.0f} dBm")
        print(f"Rx power: {self.rx_power:.0f} dBW = {db_to_dbm(self.rx_power):.0f} dBm")
        print(f"Frequency: {self.frequency / 1e6:.1f} MHz")
        print(
            f"Wavelength: {wavelength(self.frequency):.3g} m (1/4 = {wavelength(self.frequency)/4:.3g} m)"
        )
        print(f"System noise temperature: {self.temperature_system_noise:.0f} K")
        if self.temperature_antenna is not None:
            print(f" - Antenna: {self.temperature_antenna:.0f} K")
            print(f" - Line: {self.temperature_line:.0f} K")
            print(f" - Receiver: {self.temperature_receiver:.0f} K")

        g_over_t = self.rx_gain - to_db(self.temperature_system_noise)
        print(f"G/T: {g_over_t:.1f} dB/K")

        print(f"Nominal data rate: {self.data_rate/1e3:.1f} kb/s")
        if self.max_data_rate is not None:
            print(
                f"Maximum data rate (with {self.margin_target:.0f} dB margin): {self.max_data_rate/1e3:.1f} kb/s"
            )

        round_trip_time = 2 * self.latency_factor * calculate_travel_time(self.tx_rx_distance)
        print(f"Round-trip time: {round_trip_time * 1e3:.0f} ms")

    def print_table(self):
        budget = list(self.budget.itertuples(name=None))
        snr = list(self.snr.itertuples(name=None))
        print(tabulate(budget + [SEPARATING_LINE] + snr, numalign="right", floatfmt="+.3f"))
        if self.snr.loc["Margin"][0] > self.margin_target:
            print(colored("Budget closed", "green"))
        else:
            print(colored("Budget not closed", "red"))

    @classmethod
    def _get_component_value(cls, df: pd.DataFrame, component: str):
        return df.loc[component]["Value [dB]"]

    @classmethod
    def _make_plotting_df(cls):
        return pd.DataFrame(
            [
                ("Tx power", ["Tx power"], 0),
                ("Tx antenna", ["Tx line loss"], 0.15),
                ("EIRP", ["Tx antenna gain", "Tx antenna pointing loss"], 0.3),
                ("Received radiation", ["Free space loss", "Environment loss"], 0.7),
                ("Rx antenna", ["Rx antenna gain", "Rx antenna pointing loss"], 0.85),
                ("Rx power", ["Rx line loss"], 1),
            ],
            columns=["loc_name", "components", "pos"],
        )

    def plot(self, ax: plt.Axes = None):
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 4))

        df = self._make_plotting_df()

        def get_gain(row):
            gain = 0
            for component in row["components"]:
                gain += self._get_component_value(self.budget, component)
            return gain

        df["gain"] = df.apply(get_gain, axis=1)
        df["cum_gain"] = df["gain"].cumsum()

        power_required = (
            self._get_component_value(self.snr, "Required Eb/N0")
            + self.margin_target
            - self._get_component_value(self.budget, "Data rate")
            - self._get_component_value(self.budget, "Noise spectral density")
        )
        print(power_required)

        for pos in df["pos"]:
            ax.axvline(pos, color="black", alpha=0.5, ls=":")

        ax.plot(df["pos"], df["cum_gain"], c="black", label="Actual power")
        ax.hlines(
            power_required,
            0.95,
            1.05,
            ls="--",
            color="black",
            label=f"Required power (with {self.margin_target} dB margin)",
        )

    @classmethod
    def plot_multiple(cls, links: list[Self], save_name: str = None):
        _, axs = plt.subplots(
            len(links), figsize=(10, 2.5 * len(links)), sharex="all", sharey="all"
        )

        df = cls._make_plotting_df()

        for ax, link in zip(axs, links):
            link.plot(ax)
            ax.set_title(link.name)
            ax.set_ylabel("Power [dBW]")
            ax_dbm = ax.secondary_yaxis("right", functions=(db_to_dbm, dbm_to_db))
            ax_dbm.set_ylabel("Power [dBm]")

        axs[0].legend(loc="upper right")
        axs[-1].set_xticks(df["pos"], df["loc_name"])

        format_plot()
        if save_name:
            save_plot("linkbudget", save_name)
        plt.show()


if __name__ == "__main__":
    relay_gain_uhf_rx = 2
    relay_gain_uhf_tx = 5
    aircraft_gain_uhf = 4
    aircraft_gain_hf = 0
    base_gain_hf = 0

    relay_pointing_loss_uhf = -1
    aircraft_pointing_loss_uhf = -2
    aircraft_pointing_loss_hf = -3
    base_pointing_loss_hf = -3

    aircraft_power_uhf = 25
    relay_power_uhf = 15
    aircraft_power_hf = 1
    base_power_hf = 5

    tx_loss = -1
    rx_loss = -1

    base_noise_figure = 0.5
    relay_noise_figure = 4.9
    aircraft_noise_figure = 4.9

    relay_tx_rx_distance = 8240e3
    skywave_tx_rx_distance = 1050e3
    nv_skywave_tx_rx_distance = 320e3

    relay_environment_loss_uhf = -0.5
    skywave_environment_loss_hf = 2 * -0.5 - 39

    temperature_antenna_upward_looking = 54
    temperature_antenna_downward_looking = 155

    data_rate_uplink = calculate_uplink_data_rate()
    data_rate_downlink = calculate_downlink_data_rate()

    budget_uplink_relay = Link(
        name="Uplink (relay)",
        frequency=402e6,
        tx_power=aircraft_power_uhf,
        tx_loss=tx_loss,
        tx_gain=aircraft_gain_uhf,
        tx_pointing_loss=aircraft_pointing_loss_uhf,
        tx_rx_distance=relay_tx_rx_distance,
        loss_environment=relay_environment_loss_uhf,
        rx_pointing_loss=relay_pointing_loss_uhf,
        rx_gain=relay_gain_uhf_rx,
        rx_loss=rx_loss,
        noise_figure=relay_noise_figure,
        temperature_antenna=temperature_antenna_downward_looking,
        data_rate=data_rate_uplink,
        bit_error_rate=1e-6,
        modulation="bpsk",
        coding="ldpc-1/2",
        latency_factor=2,
    )
    budget_uplink_relay.print_configuration()
    budget_uplink_relay.print_table()

    print()
    budget_downlink_relay = Link(
        name="Downlink (relay)",
        frequency=437e6,
        tx_power=relay_power_uhf,
        tx_loss=tx_loss,
        tx_gain=relay_gain_uhf_tx,
        tx_pointing_loss=relay_pointing_loss_uhf,
        tx_rx_distance=relay_tx_rx_distance,
        loss_environment=relay_environment_loss_uhf,
        rx_pointing_loss=aircraft_pointing_loss_uhf,
        rx_gain=aircraft_gain_uhf,
        rx_loss=rx_loss,
        noise_figure=aircraft_noise_figure,
        temperature_antenna=temperature_antenna_upward_looking,
        data_rate=data_rate_downlink,
        bit_error_rate=1e-6,
        modulation="bpsk",
        coding="ldpc-1/2",
        latency_factor=2,
    )
    budget_downlink_relay.print_configuration()
    budget_downlink_relay.print_table()

    print()
    budget_uplink_skywave = Link(
        name="Uplink (skywave)",
        frequency=4e6,
        tx_power=aircraft_power_hf,
        tx_loss=tx_loss,
        tx_gain=aircraft_gain_hf,
        tx_pointing_loss=aircraft_pointing_loss_hf,
        tx_rx_distance=skywave_tx_rx_distance,
        loss_environment=skywave_environment_loss_hf,
        rx_pointing_loss=base_pointing_loss_hf,
        rx_gain=base_gain_hf,
        rx_loss=rx_loss,
        noise_figure=base_noise_figure,
        temperature_antenna=temperature_antenna_upward_looking,
        data_rate=data_rate_uplink,
        bit_error_rate=1e-6,
        modulation="bpsk",
        coding="ldpc-1/2",
    )
    budget_uplink_skywave.print_configuration()
    budget_uplink_skywave.print_table()

    print()
    budget_downlink_skywave = Link(
        name="Downlink (skywave)",
        frequency=4e6,
        tx_power=base_power_hf,
        tx_loss=tx_loss,
        tx_gain=base_gain_hf,
        tx_pointing_loss=base_pointing_loss_hf,
        tx_rx_distance=skywave_tx_rx_distance,
        loss_environment=skywave_environment_loss_hf,
        rx_pointing_loss=aircraft_pointing_loss_hf,
        rx_gain=aircraft_gain_hf,
        rx_loss=rx_loss,
        noise_figure=aircraft_noise_figure,
        temperature_antenna=temperature_antenna_upward_looking,
        data_rate=data_rate_downlink,
        bit_error_rate=1e-6,
        modulation="bpsk",
        coding="ldpc-1/2",
    )
    budget_downlink_skywave.print_configuration()
    budget_downlink_skywave.print_table()

    print()
    budget_uplink_nv_skywave = Link(
        name="Uplink (near-vertical skywave)",
        frequency=4e6,
        tx_power=aircraft_power_hf,
        tx_loss=tx_loss,
        tx_gain=aircraft_gain_hf,
        tx_pointing_loss=aircraft_pointing_loss_hf,
        tx_rx_distance=nv_skywave_tx_rx_distance,
        loss_environment=skywave_environment_loss_hf,
        rx_pointing_loss=base_pointing_loss_hf,
        rx_gain=base_gain_hf,
        rx_loss=rx_loss,
        noise_figure=base_noise_figure,
        temperature_antenna=temperature_antenna_upward_looking,
        data_rate=data_rate_uplink,
        bit_error_rate=1e-6,
        modulation="bpsk",
        coding="ldpc-1/2",
    )
    budget_uplink_nv_skywave.print_configuration()
    budget_uplink_nv_skywave.print_table()

    print()
    budget_downlink_nv_skywave = Link(
        name="Downlink (near-vertical skywave)",
        frequency=4e6,
        tx_power=base_power_hf,
        tx_loss=tx_loss,
        tx_gain=base_gain_hf,
        tx_pointing_loss=base_pointing_loss_hf,
        tx_rx_distance=nv_skywave_tx_rx_distance,
        loss_environment=skywave_environment_loss_hf,
        rx_pointing_loss=aircraft_pointing_loss_hf,
        rx_gain=aircraft_gain_hf,
        rx_loss=rx_loss,
        noise_figure=aircraft_noise_figure,
        temperature_antenna=temperature_antenna_upward_looking,
        data_rate=data_rate_downlink,
        bit_error_rate=1e-6,
        modulation="bpsk",
        coding="ldpc-1/2",
    )
    budget_downlink_nv_skywave.print_configuration()
    budget_downlink_nv_skywave.print_table()

    # Link.plot_multiple(
    #     [budget_uplink_relay, budget_uplink_skywave, budget_uplink_nv_skywave],
    #     "linkbudget_uplink_all",
    # )
    # Link.plot_multiple([budget_uplink_relay, budget_downlink_relay])
    # Link.plot_multiple([budget_uplink_skywave, budget_downlink_skywave])
    # Link.plot_multiple([budget_uplink_nv_skywave, budget_downlink_nv_skywave])
