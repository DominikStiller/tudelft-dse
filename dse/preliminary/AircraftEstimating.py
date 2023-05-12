import numpy as np
from tiltrotormain import constants
from cruise_sizing import area_and_thrust
from power_sizing import power
import scipy.integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from RotorEngineSizing import RadiusMassElementMomentum
#from PerformanceAnalysis import RangeCalc
from cruise_sizing import area_and_thrust
from power_sizing import power

def DragEstimation(lf, hf, Swing, t2c, Vcr, visc_cr, Cl, AR, rho):
    Oswald = 0.9

    # Fuselage
    bf = hf
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

    Cd = Cd_0 + Cl ** 2 / (np.pi * AR * Oswald)
    D = Cd * rho * 0.5 * Vcr ** 2 * Swing

    print('Drag of the aircraft without wing will be of ' + str(D) + '[N] \n')
    return D


def RangeCalc(Wto, Wtot, R, AR, V_cr, E_density, P_density, E_density_TO):
    g_mars = constants['airDensity']

    # Dimensions
    b = 2 * 1.5 * R  # Wingspan
    c = b / AR  # chord
    bf = c  # Width of the fuselage = chord
    lf = 1.78 + c * 3  # Length of the fuselage equals 1.78 + 3c
    if c < 1.78:  # Establish lower bounds
        lf = 1.78 * 4
        bf = 1.78
    hf = bf  # Height of the fuselage
    Swing = area_and_thrust(0, constants['cl'], constants['cd'], Wto, 0.5 * constants['airDensity'] * V_cr)[0]

    # Thrust estimations
    T_to = 1.1 * Wto * g_mars
    T_cr = area_and_thrust(0, constants['cl'], constants['cd'], Wto, 0.5 * constants['airDensity'] * V_cr)[1] + \
           DragEstimation(lf, hf, Swing, 0.12, V_cr, 5.167E-4, 1.2, AR, rho=0.01)

    # Power and energy

    # Take-off
    Power = power(T_to, R)
    E_TO = Power * (1 / 6)  # Assuming take-off time of 10min
    m_TO = max(Power / P_density, E_TO / E_density_TO)
    E_TO = Power * (1 / 6)  # Assuming take-off time of 10min

    # Cruise
    P_cr = power(T_cr, R)
    E_cr = ((Wto - Wtot - m_TO) * E_density)  # Use all remaining mass for batteries
    m_battery_cr = Wto - Wtot - m_TO
    Endurance = E_cr / P_cr
    Range = Endurance * V_cr * 3.6
    # Range and endurance
    Endurance = E_cr / P_cr  # hours
    Range = Endurance * V_cr * 3.6  # km

    return Range, Endurance, m_TO, m_battery_cr


# Weight Prediction:
def Class2Weight(R, RotorMass, Wto, N_ult, AR, wingbraced, V_cr, E_density, P_density, E_density_TO, m_payload,
                 m_solar):
    g_mars = constants['gravityMars']

    # Initial dimensions
    b = 2 * 1.5 * R  # wingspan
    c = b / AR  # chord
    bf = c  # Fuselage width
    hf = bf  # Fuselage height
    lf = 1.78 + c * 3  # Fuselage length
    if c < 1.78:
        lf = 1.78 * 4
        bf = 1.78

    # Wing and tail area
    Swing = area_and_thrust(0, constants['cl'], constants['cd'], Wto, 0.5 * constants['airDensity'] * V_cr)[0]
    Stail = Swing * c / (1.5 * R)

    # Wing Group
    Wwing2Wto = 4.9e-3 * b ** 0.75 * (1 + np.sqrt(1.905 / b)) * N_ult ** 0.55 * ((b / c) / (Wto / Swing)) ** 0.3
    if wingbraced:
        Wwing2Wto *= 0.7

    # Tail Group
    Wtail2Wto = 0.64 * (N_ult * Stail ** 2) ** 0.75

    # Body Group
    lamdaf = lf / hf
    Vdive = 1.25 * V_cr
    lt = lf
    S_g = np.pi * hf * lf * (1 - 2 / lamdaf) ** (2 / 3) * (1 + 1 / (lamdaf ** 2))
    Wf = .23 * np.sqrt(Vdive * lt / (bf + hf)) * S_g ** 1.2

    # Control Surfaces
    ksc = 0.64  # transport plane with powered controls
    Wsc = 0.768 * ksc * Wto ** (2 / 3)

    #Total Weight
    Wtot = Wwing2Wto*Wto+Wtail2Wto+Wf+Wsc+RotorMass+m_payload+m_solar
    Range, Endurance, m_battery_TO, m_battery_cr = RangeCalc(Wto, Wtot, R, AR, V_cr, E_density, P_density, E_density_TO)

    print('Wing weight: ' + str(Wwing2Wto * Wto) + '[kg]')
    print('Tail weight: ' + str(Wtail2Wto) + '[kg]')
    print('Body weight: ' + str(Wf) + '[kg]')
    print('Control Surfaces: ' + str(Wsc) + '[kg]')
    print('Available weight for batteries: ' + str(Wto - Wtot) + '[kg]')
    print('Available Endurance: ' + str(Endurance) + '[h]')
    print('Available Range: ' + str(Range) + '[km]')
    print('Flight Radius: ' + str(Range / 2) + '[km]')
    return Wwing2Wto * Wto, Wtail2Wto, Wf, Wsc
