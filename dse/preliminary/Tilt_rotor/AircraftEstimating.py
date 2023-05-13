from constants import const, aircraftParameters
from cruise_sizing import area_and_thrust
from power_sizing import power
import numpy as np


def DragEstimation(R, Swing, t2c, Vcr, visc_cr, AR):
    Oswald = 0.9

    # Dimensions
    b = max(np.sqrt(Swing * AR), 3*R)  # wingspan
    c = b / AR  # chord
    bf = c / 2  # Fuselage width
    lf = c * 2  # Fuselage length
    if c < 1.78:
        lf = 1.78 * 2.5
        bf = 1
    hf = bf  # Fuselage height
    # Initial dimensions

    # Fuselage
    Cd_fusS = 0.0031 * lf * (bf + hf)

    # Wing
    Cd_wingS = 0.0054 * (1 + 3 * t2c * np.cos(0) ** 2) * Swing  # Sweep is 0

    # Engines
    Cd_engineS = 2 * 0.07 * 0.08

    # Undercarriage
    r_uc = 1.08  # Main gear retracted in strealined fairings

    # Reynolds correction
    r_Re = 47 * (Vcr * lf / visc_cr) ** (-0.2)

    # Tailplane
    r_t = 1.24

    Cd_0 = r_Re * r_uc * (r_t * (Cd_fusS + Cd_engineS)) / Swing

    Cd = Cd_0 + const['cl'] ** 2 / (np.pi * AR * Oswald)
    D = Cd * const['airDensity'] * 0.5 * Vcr ** 2 * Swing

    print('Drag of the aircraft without wing will be of ' + str(D) + '[N] \n')
    return D


def RangeCalc(Wto, Wtot, R, AR, V_cr, E_density, P_density, E_density_TO):
    g_mars = const['gravityMars']

    # Dimensions

    Swing, T_wing = area_and_thrust(0, const['cl'], const['cd'], Wto, 0.5 * const['airDensity'] * V_cr ** 2)

    # Thrust estimations
    T_to = 1.1 * Wto * g_mars
    T_cr = T_wing + \
           DragEstimation(R, Swing, const['t/c'], V_cr, 5.167E-4, AR)

    # Take-off
    Power = 4 * power(T_to / 4, R)
    E_TO = Power * const['takeOffTime']/3600
    m_TO = max(Power / P_density, E_TO / E_density_TO)

    # Cruise
    P_cr = power(T_cr, R)
    E_cr = ((Wto - Wtot - m_TO) * E_density)  # Use all remaining mass for batteries
    m_battery_cr = Wto - Wtot - m_TO
    if np.min([m_battery_cr]) < 0:
        raise ValueError('Not enough available weight for batteries')

    Endurance = E_cr / P_cr
    Range = Endurance * V_cr * 3.6

    return Range, Endurance, m_TO, m_battery_cr


# Weight Prediction:
def Class2Weight(R, RotorMass, Wto, N_ult, AR, wingbraced, V_cr, E_density, P_density, E_density_TO, m_payload,
                 m_solar, print_results=True):
    Swing = area_and_thrust(0, const['cl'], const['cd'], Wto, 0.5 * const['airDensity'] * V_cr ** 2)[0]
    b = max(np.sqrt(Swing * AR), 3*R)  # wingspan
    c = b / AR  # chord
    bf = c / 2  # Fuselage width
    lf = c * 2  # Fuselage length
    if c < 1.78:
        lf = 1.78 * 2.5
        bf = 1
    hf = bf  # Fuselage height

    # Wing and tail area
    Stail = Swing * c / (1.5 * R)

    # Wing Group
    Wwing2Wto = 4.9e-3 * b ** 0.75 * (1 + np.sqrt(1.905 / b)) * N_ult ** 0.55 * ((b / c) / (Wto / Swing)) ** 0.3
    if wingbraced:
        Wwing2Wto *= 0.7

    # Tail Group
    Wtail2Wto = 0.64 * (N_ult * Stail ** 2) ** 0.75

    # Body Group
    Vdive = 1.1 * V_cr
    lt = lf
    S_g = 4/3 * np.pi * hf**3 + 2*np.pi*hf * lf
    Wf = .23 * np.sqrt(Vdive * lt / (bf + hf)) * S_g ** 1.2

    # Total Weight
    Wtot = Wwing2Wto*Wto + Wtail2Wto + Wf + RotorMass + m_payload + m_solar
    Range, Endurance, m_battery_TO, m_battery_cr = RangeCalc(Wto, Wtot, R, AR, V_cr, E_density, P_density, E_density_TO)

    if print_results:
        print('Wing weight: ' + str(Wwing2Wto * Wto) + '[kg]')
        print('Tail weight: ' + str(Wtail2Wto) + '[kg]')
        print('Body weight: ' + str(Wf) + '[kg]')
        # print('Available weight for batteries: ' + str(Wto - Wtot) + '[kg]')
        # print('Available Endurance: ' + str(Endurance) + '[h]')
        # print('Available Range: ' + str(Range) + '[km]')
        # print('Flight Radius: ' + str(Range / 2) + '[km]')
    return Range, Wwing2Wto * Wto, Wtail2Wto, Wf
