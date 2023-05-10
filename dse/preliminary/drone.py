import numpy as np
import scipy
import matplotlib.pyplot as plt

from dse.plotting import format_plot

rho = 0.01  # density
speed_of_sound = 220
M_max = 0.85  # maximum mach number
mass = 3000
g_mars = 3.71

# Geometric parameters
c_to_R_ratio = 1/20  # ratio of chord to radius
R = 17  # radius of rotor
omega = 209  # angular speed of rotor [omega 209 = 2000 RPM]
# a = 0.11 * 180 / np.pi  # lift curve slope
a = 6  # lift curve slope
theta_0 = np.radians(25)  # zero pitch angle
theta_tw = np.radians(8)  # twist angle
n_blades = 4
n_rotors = 8  # if coaxial, use number of axes
coaxial = True

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


def tip_mach(omega, R, V_fwd=0):
    V_tip = omega * R
    return (V_tip + V_fwd) / speed_of_sound


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
        chord = R * c_to_R_ratio
        SR = solidity_ratio(chord, R)
        vi = velocity_induced_hover(T, rho, A)
        V = 0
        alpha = 0
        T_calc = _thrust_single_rotor(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
        return T - T_calc

    res = scipy.optimize.least_squares(f, x0=5e3, bounds=(0, np.inf))
    T_per_rotor = res.x[0]
    if coaxial:
        return 2 * T_per_rotor * n_rotors * 0.88
    else:
        return T_per_rotor * n_rotors


def solve_thrust_forward_flight(V, alpha, R=R, omega=omega):
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
        chord = R * c_to_R_ratio
        SR = solidity_ratio(chord, R)
        vi = velocity_induced_forward_flight(T, rho, A, V)
        T_calc = _thrust_single_rotor(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
        return T - T_calc

    res = scipy.optimize.least_squares(f, x0=5e3, bounds=(0, np.inf))
    T_per_rotor = res.x[0]
    if coaxial:
        return 2 * T_per_rotor * n_rotors * 0.88
    else:
        return T_per_rotor * n_rotors


def solve_thrust_vertical_climb(V_c, R=R, omega=omega):
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
        chord = R * c_to_R_ratio
        SR = solidity_ratio(chord, R)
        vi = velocity_induced_vertical_climb(T, rho, A, V)
        T_calc = _thrust_single_rotor(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
        return T - T_calc

    res = scipy.optimize.least_squares(f, x0=5e3, bounds=(0, np.inf))
    T_per_rotor = res.x[0]
    if coaxial:
        return 2 * T_per_rotor * n_rotors * 0.88
    else:
        return T_per_rotor * n_rotors


def solve_rotor_torque_hover(R=R, omega=omega):
    V = 0
    alpha = 0
    A = area(R)
    chord = R * c_to_R_ratio
    SR = solidity_ratio(chord, R)
    T = solve_thrust_hover(R, omega)/n_rotors
    vi = velocity_induced_hover(T, rho, A)
    torquecalc = rotortorque(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
    return torquecalc


def solve_rotor_torque_forward_flight(V, alpha, R=R, omega=omega):
    """
    Calculate rotor thrust in forward flight.

    Args:
        V: forward velocity [m/s]
        alpha: angle of attack [rad]

    Returns:
        Total thrust of all rotors [N]
    """
    A = area(R)
    chord = R * c_to_R_ratio
    SR = solidity_ratio(chord, R)
    T = solve_thrust_forward_flight(V, alpha, R, omega)/n_rotors
    vi = velocity_induced_forward_flight(T, rho, A, V)
    torquecalc = rotortorque(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
    return torquecalc


def solve_rotor_torque_vertical_climb(V_c, R=R, omega=omega):
    alpha = np.pi / 2
    V = V_c
    A = area(R)
    chord = R * c_to_R_ratio
    SR = solidity_ratio(chord, R)
    T = solve_thrust_vertical_climb(V, R, omega)/n_rotors
    vi = velocity_induced_vertical_climb(T, rho, A, V)
    torquecalc = rotortorque(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
    return torquecalc

'''
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
'''

def plot_radius_rpm_range():
    RR = np.arange(0.2, 10.5, 0.1)
    rpms = np.hstack([[200, 300, 400], np.arange(500, 2501, 500)])

    R_maxs = []
    T_max = []
    omegas = []
    for rpm in rpms:
        T = []
        omega = omega_from_rpm(rpm)
        for R in RR:
            T.append(solve_thrust_hover(R, omega))
        plt.plot(RR, T, label=f"{rpm} rpm")

        R_max = M_max * speed_of_sound / omega
        R_maxs.append(R_max)
        omegas.append(omega)
        T_max.append(solve_thrust_hover(R_max, omega))

    plt.scatter(R_maxs, T_max, marker="x", label="Max. radius for M<=0.8")
    plt.xlabel("Radius [m]")
    plt.ylabel("Thrust [N]")
    plt.yscale("log")
    plt.axhline(mass * g_mars, c="black", label="Weight at 3000 kg")
    plt.legend()
    format_plot()
    plt.show()
    return R_maxs, T_max, rpms * (2 * np.pi/ 60)

def thrust_and_force_for_optimum_Tip_Mach(speed, angle_of_attack, climb_speed):
    rmax, tmaxhower, omegas = plot_radius_rpm_range()
    tmaxclimb = []
    tmaxcruise = []
    torquehower = []
    torqueclimb = []
    torquecruise = []
    for i in range(len(omegas)):
        tmaxclimb.append(solve_thrust_vertical_climb(climb_speed, rmax[i], omegas[i]))
        tmaxcruise.append(solve_thrust_forward_flight(speed, angle_of_attack, rmax[i], omegas[i]))
        torquehower.append(solve_rotor_torque_hover(rmax[i], omegas[i]))
        torqueclimb.append(solve_rotor_torque_vertical_climb(climb_speed, rmax[i], omegas[i]))
        torquecruise.append(solve_rotor_torque_forward_flight(speed, angle_of_attack, rmax[i], omegas[i]))

    return tmaxhower, tmaxclimb, tmaxcruise, torquehower, torqueclimb, torquecruise

def engine_and_fuel_mass(torque, omega):
    power = torque * omega
    hp = power/745.7
    kwh = power * 60 / 1000
    fc = (kwh * 0.025 / 0.568044) * 1.5

    # kg per hp 0.2 to 0.4 of avg helicopter engine based on PT6 engine series
    print(f"engine mass lower bound = {0.2 * hp * (1 / 0.568044)} kg per engine, {0.2 * hp * 8 * (1 / 0.568044)} kg total")
    print(f"engine mass upper bound = {0.4 * hp * (1 / 0.568044)} kg per engine, {0.4 * hp * 8 * (1 / 0.568044)} kg total")

    # fuel consumption = 205 g/kWh if 1680 hp so assume 20.5 g/kWh for 115 hp based on S6U-PTA engine
    print(f"fuel mass = {fc} kg")

if __name__ == "__main__":
    speed = 111
    angle_of_attack = np.radians(-5)
    climb_speed = 10

    rmax, tmaxhower, omegas = plot_radius_rpm_range()
    tmaxhower, tmaxclimb, tmaxcruise, torquehower, torqueclimb, torquecruise = thrust_and_force_for_optimum_Tip_Mach(speed, angle_of_attack, climb_speed)

    engine_and_fuel_mass(torquehower[0], omegas[0])
    print('thrusts')
    print(tmaxhower)
    print(tmaxclimb)
    print(tmaxcruise)
    print('torques')
    print(torquehower)
    print(torqueclimb)
    print(torquecruise)

