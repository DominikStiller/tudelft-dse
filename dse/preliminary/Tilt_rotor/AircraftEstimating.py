from .constants import const, aircraftParameters
from .power_sizing import power
import numpy as np


def DragEstimation(Swing, Vcr, visc_cr, changeDimensions=None):
    # Initial Dimensions
    lf = 1.78 + aircraftParameters['chord']  # Fuselage length
    bf = 1
    hf = bf  # Fuselage height

    if changeDimensions is not None:
        lf = changeDimensions[0]
        bf = changeDimensions[1]
        hf = bf

    # Fuselage
    Cd_fusS = 0.0031 * lf * (bf + hf)

    # Engines
    Cd_engineS = 2 * 0.07 * 0.08

    # Undercarriage
    r_uc = 1.08  # Main gear retracted in strealined fairings

    # Reynolds correction
    r_Re = 47 * (Vcr * lf / visc_cr) ** (-0.2)

    # Tailplane
    r_t = 1.24

    Cd_0 = r_Re * r_uc * (r_t * (Cd_fusS + Cd_engineS)) / Swing
    return Cd_0


def RangeCalc(Wto, Wtot, R, V_cr):
    g_mars = const['gravityMars']

    # Thrust estimations
    T_to = Wto * g_mars
    T_cr = aircraftParameters['cruiseThrust']

    # Take-off
    Power = aircraftParameters['totalRotors'] * power(T_to/aircraftParameters['totalRotors'], R)
    E_TO = Power * const['takeOffTime'] / 3600
    m_TO = max(Power / const['takeoffBatteryPowerDensity'], E_TO / const['takeoffBatteryEnergyDensity'])

    # Cruise
    # P_cr = aircraftParameters['totalRotors']*power(T_cr/aircraftParameters['totalRotors'], R)

    P_cr = aircraftParameters['cruiseThrust'] * const['cruiseSpeed']
    E_cr = (Wto - Wtot) * const['takeoffBatteryEnergyDensity'] - E_TO * 2  # Use all remaining mass for batteries

    m_battery_cr = Wto - Wtot - m_TO
    m_battery_cr = m_battery_cr*(m_battery_cr>0)
    #if np.min([m_battery_cr]) < 0:
    #    raise ValueError('Not enough available weight for batteries')

    Endurance = E_cr / P_cr
    Range = Endurance * V_cr * 3.6

    return Range, Endurance, m_TO, m_battery_cr


# Weight Prediction:
def Class2Weight(R, RotorMass, Wto, N_ult, AR, wingbraced, V_cr, E_density, P_density, E_density_TO, m_payload,
                 m_solar, print_results=True, changeDimensions=None):
    Swing = aircraftParameters['wingArea']
    b = aircraftParameters['wingspan']  # wingspan
    c = aircraftParameters['chord']  # chord
    lf = 1.78 + c
    bf = 1  # Fuselage width

    if changeDimensions is not None:
        lf = changeDimensions[0]
        bf = changeDimensions[1]

    hf = bf  # Fuselage height

    if R <= 0:
        return "Rotor radius has to be greater than zero.",0,0,0
    if RotorMass <= 0:
        return "Rotor mass has to be greater than zero.",0,0,0
    if Wto <= 0:
        return "Take-off mass has to be greater than zero.",0,0,0
    if V_cr <= 0:
        return "Cruise speed has to be greater than zero.",0,0,0

    # Wing and tail area
    Stail = Swing * c / (1.5 * R)

    # Wing Group
    Wwing2Wto = 4.9e-3 * b ** 0.75 * (1 + np.sqrt(1.905 / b)) * N_ult ** 0.55 * ((b / c) / (Wto / Swing)) ** 0.3
    if wingbraced:
        Wwing2Wto *= 0.7

    # Tail Group
    Wtail2Wto = 0.64 * (N_ult * Stail ** 2) ** 0.75

    # Body Group
    Vdive = V_cr
    lt = lf/2
    S_g = 4 * np.pi * hf**2 + 2*np.pi*hf * lf
    Wf = .23 * np.sqrt(Vdive * lt / (bf + hf)) * S_g ** 1.2

    # Total Weight
    Wtot = Wwing2Wto*Wto + Wtail2Wto + Wf + RotorMass + m_payload + m_solar
    Range, Endurance, m_battery_TO, m_battery_cr = RangeCalc(Wto, Wtot, R, V_cr)

    if print_results:
        print('Wing weight: ' + str(Wwing2Wto * Wto) + '[kg]')
        print('Tail weight: ' + str(Wtail2Wto) + '[kg]')
        print('Body weight: ' + str(Wf) + '[kg]')
        print('Available weight for batteries: ' + str(Wto - Wtot) + '[kg]')
        print('Available Endurance: ' + str(Endurance) + '[h]')
        print('Available Range: ' + str(Range) + '[km]')
        print('Flight Radius: ' + str(Range / 2) + '[km]')
    return Range, Wwing2Wto * Wto, Wtail2Wto, Wf
