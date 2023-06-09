import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy

from dse.plotting import format_plot, save_plot

rho = 0.01  # density
speed_of_sound = 220
M_max = 0.85  # maximum mach number
mass = 3000
g_mars = 3.71

# Geometric parameters
c_to_R_ratio = 1 / 20  # ratio of chord to radius
R = 7  # radius of rotor
omega = 26.2  # angular speed of rotor [omega 0.66 = 250 RPM]
# a = 0.11 * 180 / np.pi  # lift curve slope
a = 6  # lift curve slope
theta_0 = np.radians(25)  # zero pitch angle
theta_tw = np.radians(8)  # twist angle
n_blades = 6
n_rotors = 4  # if coaxial, use number of axes
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


def thrust_coefficient(T, R, omega):
    A = np.pi * R**2
    return T / (rho * A * omega**2 * R**2)


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
        V_c: climb speed [m/s]
        R: radius [rad]
        omega = angular velocity [rad/s]

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
    T = solve_thrust_hover(R, omega) / n_rotors
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
    T = solve_thrust_forward_flight(V, alpha, R, omega) / n_rotors
    vi = velocity_induced_forward_flight(T, rho, A, V)
    torquecalc = rotortorque(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
    return torquecalc


def solve_rotor_torque_vertical_climb(V_c, R=R, omega=omega):
    alpha = np.pi / 2
    V = V_c
    A = area(R)
    chord = R * c_to_R_ratio
    SR = solidity_ratio(chord, R)
    T = solve_thrust_vertical_climb(V, R, omega) / n_rotors
    vi = velocity_induced_vertical_climb(T, rho, A, V)
    torquecalc = rotortorque(rho, A, SR, a, theta_0, theta_tw, R, alpha, vi, omega, V)
    return torquecalc


"""
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
"""


def plot_radius_rpm_range():
    plt.subplots(figsize=(8, 3.5))

    RR = np.arange(0.2, 10.5, 0.1)
    rpms = np.array([200, 250, 350])

    cmap = mpl.colormaps["Blues"](np.linspace(1, 0.3, len(rpms)))

    R_maxs = []
    T_max = []
    omegas = []
    for rpm, color in zip(rpms, cmap):
        T = []
        omega = omega_from_rpm(rpm)
        for R in RR:
            T.append(solve_thrust_hover(R, omega))
        plt.plot(RR, T, label=f"{rpm} rpm", color=color)

        R_max = M_max * speed_of_sound / omega
        R_maxs.append(R_max)
        omegas.append(omega)
        T_max.append(solve_thrust_hover(R_max, omega))

    plt.scatter(R_maxs, T_max, marker="x", c="black", label="Max. radius for $M \leq 0.85$")
    plt.axhline(mass * g_mars, c="darkslategray", ls=":", label="Weight at 3000 kg")

    plt.xlabel("Blade radius [m]")
    plt.ylabel("Total thrust [N]")
    plt.yscale("log")
    plt.xlim([0.5, plt.xlim()[1]])
    plt.ylim([10, plt.ylim()[1]])

    plt.legend(loc="lower right", ncol=2)

    format_plot(ygrid=False)
    save_plot(".", "multicopter_radius_rpm")

    plt.show()
    return R_maxs, T_max, rpms * (2 * np.pi / 60)


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
        torquecruise.append(
            solve_rotor_torque_forward_flight(speed, angle_of_attack, rmax[i], omegas[i])
        )

    return tmaxhower, tmaxclimb, tmaxcruise, torquehower, torqueclimb, torquecruise


def engine_and_fuel_mass(torque, omega):
    power = torque * omega
    hp = power / 745.7
    kwh = power * 60 / 1000
    fc = (kwh * 0.025 / 0.568044) * 2  # *2 added as safety factor
    # kg per hp 0.2 to 0.4 of avg helicopter engine based on PT6 engine series
    print(
        f"engine mass lower bound = {0.2 * hp * (1 / 0.568044)} kg per engine, {0.2 * hp * 4 * (1 / 0.568044)} kg total"
    )
    print(
        f"engine mass upper bound = {0.4 * hp * (1 / 0.568044)} kg per engine, {0.4 * hp * 4 * (1 / 0.568044)} kg total"
    )

    # fuel consumption = 205 g/kWh if 1680 hp so assume 20.5 g/kWh for 115 hp based on S6U-PTA engine
    print(f"fuel mass = {fc} kg")


def dragandliftofbody(V):
    q = 0.5 * rho * V**2

    # fus = 3 * 2 * 5 area approx
    aoafus = [-8, -6, -4, -2, 0, 2, 4, 6, 8]
    CDfus = [0.10, 0.09, 0.10, 0.10, 0.10, 0.11, 0.11, 0.12, 0.13]
    CLfus = [-0.06, -0.04, -0.02, 0, 0.02, 0.03, 0.04, 0.06, 0.08]
    Afus = 3 * 2

    # blade drag coefficients
    alphablade = aoafus[2] + 8 + 17 * (1 / 2)
    # acording to koning 2019 at M = 0.5
    clablade = 0.1
    cl_cd = 15
    Clblade = clablade * alphablade
    Cdblade = Clblade / cl_cd
    lblades = 7
    ARb = 20
    e = 0.7
    nrblades = 4
    nrrotors = 4  # not rotors but like nr of coaxial rotors

    # connector drag coefficient
    Cdconnect = 0.04
    Clconnect = 0.6
    lconnect = 7 * 1.2
    ARc = 20

    Dfus = (CDfus[2] + (CLfus[2] ** 2) / (np.pi * (5 / 3) * e)) * Afus * q
    Dblades = (
        (Cdblade + (Clblade**2) / (np.pi * ARb * e))
        * nrblades
        * nrrotors
        * ((lblades**2) / ARb)
        * q
    )
    Dconnect = (
        (Cdconnect + (Clconnect**2) / (np.pi * ARc * e)) * nrrotors * ((lconnect**2 / ARc)) * q
    )

    Lfus = CLfus[2] * q * Afus
    Lconnect = Clconnect * nrrotors * ((lconnect**2 / ARc)) * q

    print(
        f"The total drag is: {Dfus + Dblades + Dconnect} \n The drag of the blades is {Dblades} \n The drag of the connectors is {Dconnect} \n The drag of the fuselage is {Dfus}"
    )
    print(
        f"The total lift is: {Lfus + Lconnect} \n The lift of the connectors is {Lconnect} \n The lift of the fuselage is {Lfus}"
    )


if __name__ == "__main__":
    ct = thrust_coefficient(solve_thrust_hover(R, omega) / (2 * n_rotors), R, omega)
    sigma = solidity_ratio(R * c_to_R_ratio, R)
    print(ct / sigma, sigma)

    print(solve_thrust_hover(R, 26))
    alpha = 5
    print(solve_thrust_forward_flight(111, np.radians(alpha), R, 11) * np.cos(np.radians(alpha)))
    print(solve_thrust_forward_flight(111, np.radians(alpha), R, 11) * np.sin(np.radians(alpha)))

    speed = 111
    angle_of_attack = np.radians(10)
    climb_speed = 10

    plot_radius_rpm_range()

    # rmax, tmaxhower, omegas = plot_radius_rpm_range()
    # tmaxhower, tmaxclimb, tmaxcruise, torquehower, torqueclimb, torquecruise = thrust_and_force_for_optimum_Tip_Mach(speed, angle_of_attack, climb_speed)
    #
    # engine_and_fuel_mass(torquehower[0], omegas[0])
    #
    # dragandliftofbody(speed)
    #
    # print('thrusts')
    # print(tmaxhower)
    # print(tmaxclimb)
    # print(tmaxcruise)
    # print('torques')
    # print(torquehower)
    # print(torqueclimb)
    # print(torquecruise)
