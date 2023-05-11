import numpy as np
import matplotlib.pyplot as plt


def area_and_thrust(thrust_deflection, cl, cd, MTOM, q):
    g_mars=3.721
    A = np.array([[cl * q, np.sin(np.radians(thrust_deflection))],
                  [-cd * q, np.cos(np.radians(thrust_deflection))]])
    B = np.array([MTOM * g_mars, 0])
    S, T = np.linalg.solve(A, B)
    return S, T


if __name__ == '__main__':
    MTOM = 3000
    g_mars = 3.71
    max_thrust = 1.1 * MTOM * g_mars

    # Assumptions
    # 2D wing
    # Fly at AoA for Cl_max
    # Thrust is constant


    # OPTn airfoil data - AOA = 8 deg
    cl = 1.5
    cd = 0.11

    # Obtain data for S and gamma_e
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all')
    fig.subplots_adjust(hspace=0)

    rho = 0.01
    for V in [112, 0.7*220]:
        q = 0.5*rho*V**2
        results_list = list()
        for y in range(0, 91):
            S, T = area_and_thrust(y, cl, cd, MTOM)
            results_list.append([y, S, T])

        results_arr = np.array(results_list).T

        ax1.plot(results_arr[0], results_arr[1], label=f'V = {round(V, 1)}')
        ax2.plot(results_arr[0], results_arr[2])

    ax1.set(ylabel=r'Wing area [m$^2$]')
    ax1.legend()
    ax1.grid()
    ax2.plot([results_arr[0][0], results_arr[0][-1]], [max_thrust, max_thrust])
    ax2.set(xlabel='Thrust deflection [deg]', ylabel='Thrust [N]')
    ax2.grid()
    #plt.savefig('S&T-vs-gamma.png')
    plt.show()
