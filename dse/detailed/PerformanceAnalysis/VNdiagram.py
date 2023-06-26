import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dse.plotting import format_plot
from dse.plotting import save_plot


def n_load_quadratic(V):
    rho = 0.01
    C_L_max = 2.4
    m = 3000
    g_M = 3.71
    W = m * g_M
    S = 133.5
    n_load_quadratic = (0.5 * rho * V**2 * C_L_max) / (W / S)
    return n_load_quadratic


if __name__ == "__main__":
    # V-n diagram
    V = np.linspace(0, 138.75, 1000)
    n = n_load_quadratic(V)
    n[np.where(n > 2.2)] = 2.2
    n2 = -1 * n_load_quadratic(V)
    n2[np.where(n2 < -1)] = -1
    n2[np.where(V > 111.11)] = np.linspace(-1, 0, int(np.size(V[np.where(V > 111.11)[0][0] :])))

    # Gust diagram
    n_gust = [1, 1.5217, 1.531, 1.329, 0.671, 0.469, 0.478, 1]
    V_gust = [0, 83.35, 112, 138.75, 138.75, 112, 83.35, 0]

    plt.plot(V, n, color="b", label="Load envelope")
    plt.plot(V, n2, color="b")
    plt.vlines(138.75, 0, 2.2, color="b")
    plt.plot(V_gust, n_gust, color="darkgray", linestyle="--", label="Gust envelope")
    plt.axvline(111.11, color="black", linestyle="-.", label=r"V$_{cr}$")
    plt.ylabel("Load factor [-]")
    plt.xlabel("Velocity [m/s]")
    plt.legend()
    format_plot()
    save_plot(".", "Vn_diagram")
    plt.show()
