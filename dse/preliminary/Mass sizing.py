import numpy as np
import matplotlib.pyplot as plt
from RotorEngineSizing import RadiusfromMass


# Function calculates weight of the blades, engines, landing gear, fuselage and propulsion system
# Function does not account for payload, fuel, hydraulics, instruments, electrical system, avionics or cockpit controls
def weights(n_blades, n_legs, n_engines, radius, tip_speed, MTOM, wet_area, engine_mass):
    # Convert to imperial
    c, R, OR = np.array([radius/25, radius, tip_speed]) * 3.281
    GW, We = np.array([MTOM, engine_mass]) * 2.205
    Sw = wet_area * 3.281**2
    Lf = 1.5 * R

    g_E = 9.81 * 3.281
    J = R/25 * R**3 / 3

    # Obtain weights in lbs
    W_main_rotor_blades = 0.026 * n_blades**(2/3) * c * R**1.3 * OR**(2/3)
    W_main_rotor_hub_and_hinge = 0.0037 * n_blades**0.28 * R**1.5 * OR**0.43 * (2/3 * W_main_rotor_blades + g_E*J/R**2)**0.55
    W_fuselage = 6.9 * (GW / 1000)**0.49 * Lf**0.61 * Sw**0.25
    W_legs = 1.1 * 40 * (GW / 1000)**(2/3) * n_legs**0.54
    W_engines = engine_mass * n_engines
    W_misc = 200 * 2.205

    return np.array([4 * W_main_rotor_blades, 2 * W_main_rotor_hub_and_hinge, W_fuselage, W_legs, W_engines, W_misc]) / 2.205


def radius_mass_iteration():
    diff1 = 10
    useful_mass = 800  # Payload + fuel
    tip_speed = 200
    m_e = 122  # Dummy
    MTOM = 2000

    m_list = list()
    r_list = list()
    while diff1 > 0.01:
        R = RadiusfromMass(MTOM)[0]
        r_list.append(R)

        diff2 = 10
        Sw = 2 * np.pi * 1 * 1.5 * R + np.pi * 1 ** 2  # Cylinder of radius 1m
        while diff2 > 0.01:
            W_array = weights(n_blades=6, n_engines=2, n_legs=2,
                              radius=R, tip_speed=tip_speed, MTOM=MTOM,
                              wet_area=Sw, engine_mass=m_e)
            diff2 = abs(((np.sum(W_array) + useful_mass) - MTOM) / MTOM)
            MTOM = np.sum(W_array) + useful_mass

        m_list.append(MTOM)

        diff1 = (RadiusfromMass(MTOM)[0] - R) / R

        if len(m_list) > 8:
            plt.scatter(np.array(m_list) / 1000, r_list)
            plt.xlabel('MTOM [T]')
            plt.ylabel('R [m]')
            plt.show()
            return OverflowError

    return MTOM, R


def aircraft_estimation(Rf, Lf, Lt, MTOM, LD, Wpress):
    n_ult = 4.4
    g_M = 3.71
    rho = 0.01
    v_cruise = 154

    q = 0.5 * rho * v_cruise**2
    Sf = 2*np.pi*Rf * Lf
    fuselage_weight = 0.052 * Sf**1.086 * (n_ult * MTOM * g_M)**0.177 * Lt**(-0.051) * LD**(-0.072) * q**0.241 + Wpress
    return fuselage_weight


print(aircraft_estimation(1, 15, 13, 3000, 1.5/0.11, 0))
