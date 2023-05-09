import matplotlib.pyplot as plt
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
chord = 0.2
a = 0.1 * 180 / np.pi  # lift curve slope
theta_0 = np.radians(15)  # zero pitch angle
theta_tw = np.radians(5)  # twist angle
n_blades = 3
n_rotors = 4

Cd = 0.03  # avg drag coefficient


def area(R):
    return np.pi * R**2


def solidity_ratio(c, R):
    return n_blades * c / (np.pi * R)


def _thrust_single_blade(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V):
    t1 = (1 / 6) * rho * A * SR * a * theta_0 * R**2 + (1 / 8) * SR * a * rho * A * theta_tw * R**2
    t2 = (1 / 4) * rho * A * SR * a * theta_0 + (1 / 8) * SR * a * rho * A * theta_tw
    t3 = (-1 / 4) * SR * a * rho * A * R
    T = (omega**2) * t1 + ((V * np.cos(alpha)) ** 2) * t2 + omega * (vi + V * np.sin(alpha)) * t3
    return T

def rotortorque(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V):
    q1 = (1 / 8) * rho * A * SR * Cd * R**3
    q2 = (1 / 8) * rho * A * SR * Cd
    q3 = (1 / 4) * SR * a * rho * A * theta_0 * R**2 + (1 / 8) * rho * A * SR * a * theta_tw * R**2
    q4 = (-1 / 4) * SR * a * rho * A * R
    Q = (
        (omega**2) * q1
        + ((V * np.cos(alpha)) ** 2) * q2
        + omega * (vi + V * np.sin(alpha)) * q3
        + ((vi + V * np.sin(alpha)) ** 2) * q4
    )
    return Q

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

def solve_thrust_vertical_climb(V_c):
    """
    Calculate rotor thrust in forward flight.

    Args:
        V: forward velocity [m/s]
        alpha: angle of attack [rad]

    Returns:
        Thrust [N]
    """
    alpha = np.pi/2
    V = 0 #V_c
    def f(T):
        A = area(R)
        SR = solidity_ratio(chord, R)
        vi = velocity_induced_vertical_climb(T, rho, A, V)
        T_calc = _thrust_single_blade(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V_c)
        print('---------------')
        print(T)
        print(T_calc)
        return T - T_calc

    res = scipy.optimize.least_squares(f, x0=12, bounds=(-1000, np.inf))
    T_per_blade = res.x[0]
    return T_per_blade * n_blades * n_rotors

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
            Thrust [N]
        """
    A = area(R)
    SR = solidity_ratio(chord, R)
    T = solve_thrust_forward_flight(V, alpha)
    vi = velocity_induced_forward_flight(T, rho, A, V)
    torquecalc = rotortorque(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
    return torquecalc

def solve_rotor_torque_vertical_climb(V_c):
    """
        Calculate rotor thrust in forward flight.

        Args:
            V: forward velocity [m/s]
            alpha: angle of attack [rad]

        Returns:
            Thrust [N]
        """
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


if __name__ == "__main__":
    speed = 111
    angle_of_attack = -0.05
    climb_speed = 100
    print("forces:")
    #print(solve_thrust_hover())
    #print(solve_thrust_forward_flight(speed, angle_of_attack))
    print(solve_thrust_vertical_climb(climb_speed))
    print("torques")
    #print(solve_rotor_torque_hover())
    #print(solve_rotor_torque_forward_flight(speed, angle_of_attack))
    #print(solve_rotor_torque_vertical_climb(climb_speed))



'''
    alpha = np.arange(-60, 60, 0.1)
    thrust = []
    for i in alpha:
        thrust.append(solve_thrust_forward_flight(111, np.radians(i)))

    plt.plot(alpha, thrust)
    plt.show()
'''