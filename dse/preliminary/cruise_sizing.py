import numpy as np
import scipy.stats as sstats
import matplotlib.pyplot as plt


def area_and_thrust(thrust_deflection, cl, cd, MTOM, q):
    g_mars=3.721
    A = np.array([[cl * q, np.sin(np.radians(thrust_deflection))],
                  [-cd * q, np.cos(np.radians(thrust_deflection))]])
    B = np.array([MTOM * g_mars, 0])
    S, T = np.linalg.solve(A, B)
    return S, T


def max_tipSpeed(cruiseVelocity):
    maxTipSpeed = np.sqrt((0.92 * 220) ** 2 - cruiseVelocity ** 2)
    tipSpeedVsCruiseSpeed = sstats.linregress(cruiseVelocity, maxTipSpeed)
    return tipSpeedVsCruiseSpeed.slope, tipSpeedVsCruiseSpeed.intercept


if __name__ == '__main__':
    MTOM = 3000
    g_mars = 3.71
    max_thrust = 1.1 * MTOM * g_mars

    # Assumptions
    # 2D wing
    # Fly at AoA for Cl_max
    # Thrust is constant


    # Airfoil data
    cl = 1.6
    cd = 0.05

    # Obtain data for S and gamma_e
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all')
    fig.subplots_adjust(hspace=0)

    rho = 0.01
    velocityArray = np.linspace(112, 154, 25)
    for i in range(len(velocityArray)):
        V = velocityArray[i]
        q = 0.5*rho*V**2
        results_list = list()
        for y in range(0, 91):
            S, T = area_and_thrust(y, cl, cd, MTOM, q)
            results_list.append([y, S, T])

        results_arr = np.array(results_list).T

        if i % 5 == 0:
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

    # Calculate tip velocity
    a, b = max_tipSpeed(velocityArray)
    print(f'Tip speed = {np.round(a, 2)}cruiseSpeed + {np.round(b, 2)}')

