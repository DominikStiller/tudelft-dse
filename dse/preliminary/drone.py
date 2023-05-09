import numpy as np
import scipy


def velocity_induced_hover(T, rho, A):
    return np.sqrt(T / (2 * rho * A))


def velocity_induced_vertical_climb(T, rho, A, V_c):
    return -V_c / 2 + np.sqrt((V_c / 2) ** 2 + T / (2 * rho * A))


def velocity_induced_forward_flight(T, rho, A, V):
    return np.sqrt(-(V**2) / 2 + np.sqrt(V**4 / 4 + (T / (2 * rho * A)) ** 2))


rho = 0.01  # density
R = 10  # radius of rotor
omega = 4 * (2 * np.pi)  # anglular speed of rotor
# a = 2 * np.pi               # lift curve slope
chord = 0.4
a = 0.1 * 180 / np.pi  # lift curve slope
theta_0 = np.radians(8)  # zero pitch angle
theta_tw = np.radians(8)  # twist angle
n_blades = 3
n_rotors = 4

Cd = 0.01  # avg drag coefficient


def area(R):
    return np.pi * R**2


def solidity_ratio(c, R):
    return n_blades * c / (np.pi * R)


def _thrust_single_blade(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V):
    t1 = (1 / 6) * rho * A * SR * a * theta_0 * R**2 + (
        1 / 8
    ) * SR * a * rho * A * theta_tw * R**2
    t2 = (1 / 4) * rho * A * SR * a * theta_0 + (1 / 8) * SR * a * rho * A * theta_tw
    t3 = (-1 / 4) * SR * a * rho * A * R
    T = (omega**2) * t1 + (V * np.cos(alpha)) ** 2 * t2 + omega * (vi + V * np.sin(alpha)) * t3
    return T


def solve_thrust_hover():
    def f(T):
        A = area(R)
        SR = solidity_ratio(chord, R)
        vi = velocity_induced_hover(T, rho, A)
        V = 0
        alpha = 0
        T_calc = _thrust_single_blade(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
        return T - T_calc

    res = scipy.optimize.least_squares(f, x0=5e3, bounds=(1, np.inf))
    T_per_blade = res.x[0]
    return T_per_blade * n_blades * n_rotors


def solve_thrust_forward_flight(V, alpha):
    """
    Calculate rotor thrust in forward flight.

    Args:
        V: forward velocity [m/s]
        alpha: angle of attack [rad]

    Returns:
        Thrust [N]
    """

    def f(T):
        A = area(R)
        SR = solidity_ratio(chord, R)
        vi = velocity_induced_forward_flight(T, rho, A, V)
        T_calc = _thrust_single_blade(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
        return T - T_calc

    res = scipy.optimize.least_squares(f, x0=5e3, bounds=(1, np.inf))
    T_per_blade = res.x[0]
    return T_per_blade * n_blades * n_rotors


def hubforce():
    h1 = (1 / 4) * rho * A * SR * Cd * R
    h2 = (1 / 4) * rho * A * SR * a * (theta_0 + theta_tw / 2)
    H = (omega * V * np.cos(alpha)) * h1 + V * np.cos(alpha) * (vi + V * np.sin(alpha)) * h2
    return H


def rotortorque():
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


def rollingmoment():
    r1 = (-1 / 6) * SR * a * rho * A * theta_0 * R**2 - (
        1 / 8
    ) * rho * A * SR * a * theta_tw * R**2
    r2 = (1 / 8) * SR * a * rho * A * R
    R = (omega * V * np.cos(alpha)) * r1 + (V * np.cos(alpha))(vi + V * np.sin(alpha)) * r2
    return R


if __name__ == "__main__":
    print(solve_thrust_hover())
    print(solve_thrust_forward_flight(111, 0))
