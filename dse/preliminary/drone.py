import numpy as np
import scipy
import matplotlib.pyplot as plt

from dse.plotting import format_plot

rho = 0.01  # density
speed_of_sound = 220
M_max = 0.8  # maximum mach number

# Geometric parameters
R = 3  # radius of rotor
omega = 209  # angular speed of rotor [omega 209 = 2000 RPM]
chord = 0.3
a = 0.1 * 180 / np.pi  # lift curve slope
theta_0 = np.radians(15)  # zero pitch angle
theta_tw = np.radians(5)  # twist angle
n_blades = 5
n_rotors = 4

Cd = 0.03  # avg drag coefficient


def velocity_induced_hover(T_per_rotor, rho, A):
    return np.sqrt(T_per_rotor / (2 * rho * A))


def velocity_induced_vertical_climb(T, rho, A, V_c):
    return -V_c / 2 + np.sqrt((V_c / 2) ** 2 + T / (2 * rho * A))


def velocity_induced_forward_flight(T, rho, A, V):
    return np.sqrt(-(V**2) / 2 + np.sqrt(V**4 / 4 + (T / (2 * rho * A)) ** 2))


def area(R):
    return np.pi * R**2


def solidity_ratio(c, R):
    return n_blades * c / (np.pi * R)


def tip_mach(omega, R):
    v_tip = omega * R
    return v_tip / speed_of_sound


def omega_from_rpm(rpm):
    return rpm / 60 * 2 * np.pi


def _thrust_single_rotor(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V):
    t1 = (1 / 6) * rho * A * SR * a * theta_0 * R**2 + (
        1 / 8
    ) * SR * a * rho * A * theta_tw * R**2
    t2 = (1 / 4) * rho * A * SR * a * theta_0 + (1 / 8) * SR * a * rho * A * theta_tw
    t3 = (-1 / 4) * SR * a * rho * A * R
    T = (omega**2) * t1 + ((V * np.cos(alpha)) ** 2) * t2 + omega * (vi + V * np.sin(alpha)) * t3
    return T


def rotortorque(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V):
    q1 = (1 / 8) * rho * A * SR * Cd * R**3
    q2 = (1 / 8) * rho * A * SR * Cd
    q3 = (1 / 4) * SR * a * rho * A * theta_0 * R**2 + (
        1 / 8
    ) * rho * A * SR * a * theta_tw * R**2
    q4 = (-1 / 4) * SR * a * rho * A * R
    Q = (
        (omega**2) * q1
        + ((V * np.cos(alpha)) ** 2) * q2
        + omega * (vi + V * np.sin(alpha)) * q3
        + ((vi + V * np.sin(alpha)) ** 2) * q4
    )
    return Q


def solve_thrust_hover(R=R, omega=omega):
    def f(T):
        A = area(R)
        SR = solidity_ratio(chord, R)
        vi = velocity_induced_hover(T, rho, A)
        V = 0
        alpha = 0
        T_calc = _thrust_single_rotor(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
        return T - T_calc

    res = scipy.optimize.least_squares(f, x0=5e3, bounds=(0, np.inf))
    T_per_rotor = res.x[0]
    return T_per_rotor * n_rotors


def solve_thrust_forward_flight(V, alpha):
    """
    Calculate rotor thrust in forward flight.

    Args:
        V: forward velocity [m/s]
        alpha: angle of attack [rad]

    Returns:
        Total thrust of all rotors [N]
    """

    def f(T):
        A = area(R)
        SR = solidity_ratio(chord, R)
        vi = velocity_induced_forward_flight(T, rho, A, V)
        T_calc = _thrust_single_rotor(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
        return T - T_calc

    res = scipy.optimize.least_squares(f, x0=5e3, bounds=(0, np.inf))
    T_per_rotor = res.x[0]
    return T_per_rotor * n_rotors


def solve_thrust_vertical_climb(V_c):
    """
    Calculate rotor thrust in forward flight.

    Args:
        V: forward velocity [m/s]
        alpha: angle of attack [rad]

    Returns:
        Total thrust of all rotors [N]
    """
    alpha = np.pi / 2
    V = V_c

    def f(T):
        A = area(R)
        SR = solidity_ratio(chord, R)
        vi = velocity_induced_vertical_climb(T, rho, A, V)
        T_calc = _thrust_single_rotor(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
        return T - T_calc

    res = scipy.optimize.least_squares(f, x0=5e3, bounds=(0, np.inf))
    T_per_rotor = res.x[0]
    return T_per_rotor * n_rotors


def solve_rotor_torque_hover():
    V = 0
    alpha = 0
    A = area(R)
    SR = solidity_ratio(chord, R)
    T = solve_thrust_hover()
    vi = velocity_induced_hover(T, rho, A)
    torquecalc = rotortorque(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
    return torquecalc


def solve_rotor_torque_forward_flight(V, alpha):
    """
    Calculate rotor thrust in forward flight.

    Args:
        V: forward velocity [m/s]
        alpha: angle of attack [rad]

    Returns:
        Total thrust of all rotors [N]
    """
    A = area(R)
    SR = solidity_ratio(chord, R)
    T = solve_thrust_forward_flight(V, alpha)
    vi = velocity_induced_forward_flight(T, rho, A, V)
    torquecalc = rotortorque(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
    return torquecalc


def solve_rotor_torque_vertical_climb(V_c):
    alpha = np.pi / 2
    V = V_c
    A = area(R)
    SR = solidity_ratio(chord, R)
    T = solve_thrust_vertical_climb(V)
    vi = velocity_induced_vertical_climb(T, rho, A, V)
    torquecalc = rotortorque(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
    return torquecalc


def hubforce():
    h1 = (1 / 4) * rho * A * SR * Cd * R
    h2 = (1 / 4) * rho * A * SR * a * (theta_0 + theta_tw / 2)
    H = (omega * V * np.cos(alpha)) * h1 + V * np.cos(alpha) * (vi + V * np.sin(alpha)) * h2
    return H


def rollingmoment():
    r1 = (-1 / 6) * SR * a * rho * A * theta_0 * R**2 - (
        1 / 8
    ) * rho * A * SR * a * theta_tw * R**2
    r2 = (1 / 8) * SR * a * rho * A * R
    R = (omega * V * np.cos(alpha)) * r1 + (V * np.cos(alpha))(vi + V * np.sin(alpha)) * r2
    return R


def plot_radius_rpm_range():
    RR = np.arange(0.2, 5.5, 0.1)
    rpms = np.hstack([[200], np.arange(500, 2501, 500)])

    R_maxs = []
    T_max = []
    for rpm in rpms:
        T = []
        omega = omega_from_rpm(rpm)
        for R in RR:
            T.append(solve_thrust_hover(R, omega))
        plt.plot(RR, T, label=f"{rpm} rpm")

        R_max = M_max * speed_of_sound / omega
        R_maxs.append(R_max)
        T_max.append(solve_thrust_hover(R_max, omega))

    plt.scatter(R_maxs, T_max, marker="x", label="Max. radius for M<=0.8")
    plt.xlabel("Radius [m]")
    plt.ylabel("Thrust [N]")
    plt.yscale("log")
    plt.axhline(3.71 * 3000, c="black", label="Weight at 3000 kg")
    plt.legend()
    format_plot()
    plt.show()


if __name__ == "__main__":
    speed = 111
    angle_of_attack = -0.05
    climb_speed = 100
    print("forces:")
    # print(solve_thrust_hover())
    # print(solve_thrust_forward_flight(speed, angle_of_attack))
    # print(solve_thrust_vertical_climb(20))
    print("torques")
    # print(solve_rotor_torque_hover())
    # print(solve_rotor_torque_forward_flight(speed, angle_of_attack))
    # print(solve_rotor_torque_vertical_climb(climb_speed))

    plot_radius_rpm_range()


"""
    alpha = np.arange(-60, 60, 0.1)
    thrust = []
    for i in alpha:
        thrust.append(solve_thrust_forward_flight(111, np.radians(i)))

    plt.plot(alpha, thrust)
    plt.show()
"""
