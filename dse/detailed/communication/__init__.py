from math import log10

k = 1.380649e-23  # J/K
c = 299792458  # m/s


def to_db(x: float):
    return 10 * log10(x)


def from_db(x: float):
    return 10 ** (x / 10)
