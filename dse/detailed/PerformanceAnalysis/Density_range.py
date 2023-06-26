import matplotlib.pyplot as plt
import numpy as np
from plotting import format_plot, save_plot
import scipy.interpolate as sc


if __name__ == "__main__":
    TO_energy = 159206.51
    energy_available = 924.875 * 1e3 - 2 * TO_energy
    power_required = 2 * (
        np.array([15723.13174, 19653.91467, 23584.69761, 27515.48054, 31446.26347])
        + np.array([91515.96522, 114394.9565, 110332.6575, 93312.67935, 106643.0621])
    )
    max_range = 400 * (energy_available / power_required)
    density = np.linspace(0.01, 0.02, np.size(power_required))

    d2 = np.linspace(0.01, 0.02, 25)
    r2 = 4e6 * d2**2 - 144043 * d2 + 2175.6

    rel = sc.PchipInterpolator(
        np.array([0.01, 0.015, 0.02]), np.array([max_range[0], max_range[2], max_range[-1]])
    )
    r3 = rel(d2)

    plt.figure(figsize=(9, 5))
    # plt.plot(density, max_range.T, label='Data points')
    # plt.plot(d2, r2, label='Least-squares')
    plt.plot(d2, r3, label="3 points")
    # plt.plot(d3, r4, label='tweaked')

    plt.xlabel(r"Atmospheric density [kg/m$^3$]")
    plt.ylabel("Range [km]")
    # plt.legend()
    format_plot()
    save_plot(".", "density_range_diagram")
    plt.show()

    h = np.linspace(0, 5000, 100)
    p = 0.699 * np.exp(-0.00009 * h)
    T = -31 - 0.000998 * h
    rho = p / (0.1921 * (T + 273.15))

    plt.figure(figsize=(9, 5))
    plt.plot(h / 1e3, rho)
    plt.hlines(0.011, xmin=0, xmax=5, colors="black", linestyles="--")
    plt.xlabel("Altitude [km]")
    plt.ylabel(r"Density [kg/m$^3$]")
    format_plot()
    save_plot(".", "altitude_density_diagram")
    plt.show()
