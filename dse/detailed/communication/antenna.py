from math import pi, sqrt

from dse.detailed.communication.waves import wavelength


def calculate_parabolic_gain(diameter: float, frequency: float, efficiency: float):
    # diameter [m]
    # frequency [Hz]
    # efficiency [-]
    return (pi * diameter / wavelength(frequency)) ** 2 * efficiency


def calculate_parabolic_half_power_angle(diameter: float, frequency: float):
    # diameter [m]
    # frequency [Hz]
    return 21 / (frequency / 1e9 * diameter)


def calculate_helical_gain(diameter: float, length: float, frequency: float, efficiency: float):
    # diameter [m]
    # length [m]
    # frequency [Hz]
    # efficiency [-]
    return (pi * diameter) ** 2 * length / wavelength(frequency) ** 3 * efficiency


def calculate_helical_half_power_angle(diameter: float, length: float, frequency: float):
    # diameter [m]
    # length [m]
    # frequency [Hz]
    return 52 / sqrt(pi**2 * diameter**2 * length / wavelength(frequency) ** 3)


def pointing_loss(pointing_offset: float, half_power_angle: float):
    # diameter [m]
    # frequency [GHz]
    # pointing_offset [deg]
    # half_power_angle [deg]
    return -12 * (pointing_offset / half_power_angle) ** 2
