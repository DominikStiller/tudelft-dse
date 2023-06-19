from math import pi

from dse.detailed.communication import c


def wavelength(frequency: float):
    # frequency [Hz]
    # wavelength [m]
    return c / frequency


def calculate_free_space_loss(distance: float, frequency: float):
    # distance [m]
    # frequency [Hz]
    return (wavelength(frequency) / (4 * pi * distance)) ** 2


def calculate_travel_time(distance: float):
    # distance [m]
    return distance / c
