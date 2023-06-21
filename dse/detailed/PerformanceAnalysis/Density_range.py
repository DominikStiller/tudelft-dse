import matplotlib.pyplot as plt
import numpy as np
from plotting import format_plot, save_plot


if __name__ == '__main__':
    TO_energy = 59206.51
    energy_available = np.array([[790.125], [924.875]]) * 1e3 - 2 * TO_energy
    power_required = 4 * np.ones((2, 1)) * np.array([46589.11257, 103365.1217, 113701.6338, 124038.146, 134374.6582,
                                                 144711.1703, 124648.1978, 132958.0776, 141267.9575, 149577.8373])
    max_range = 400 * (energy_available / power_required)
    density = np.linspace(0.01, 0.02, np.shape(power_required)[1])

    plt.figure(figsize=(9, 2.5))
    plt.plot(density, max_range.T)
    plt.legend(['Helios', 'Atlas'])
    plt.xlabel(r'Atmospheric density [kg/m$^3$]')
    plt.ylabel('Range [km]')
    format_plot()
    save_plot('.', 'density_range_diagram')
    plt.show()
