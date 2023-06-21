import matplotlib.pyplot as plt
import numpy as np
from plotting import format_plot, save_plot


if __name__ == '__main__':
    energy_available = np.array([[790.125], [924.875]]) * 1e3
    flight_time = 1000/400  # hours
    power_required = ...
    max_range = energy_available / (flight_time * power_required)
    density = np.arange(0.01, 0.02, np.shape(power_required)[1])

    plt.figure(figsize=(9, 2.5))
    plt.plot(density, max_range)
    plt.legend(['Helios', 'Atlas'])
    plt.xlabel(r'Atmospheric density [kg/m$^3$]')
    plt.ylabel('Range [km]')
    format_plot()
    save_plot('.', 'density_range_diagram')
    plt.show()
