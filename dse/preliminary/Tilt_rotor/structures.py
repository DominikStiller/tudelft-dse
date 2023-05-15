from constants import const, aircraftParameters
import numpy as np


def calculate_cg():
    cockpitLength = 1.78
    fuselageLength = cockpitLength * 2.5

    masses = np.array([aircraftParameters['rotorMass'],
                       aircraftParameters['wingMass'],
                       aircraftParameters['bodyMass'],
                       aircraftParameters['tailMass'],
                       aircraftParameters['batteryMass'],
                       aircraftParameters['panelMass'],
                       const['gravityMars'] * 400
                       ])

    locations = np.array([cockpitLength + aircraftParameters['chord'] / 2,
                          cockpitLength + aircraftParameters['chord'] / 2,
                          fuselageLength / 2,
                          cockpitLength + aircraftParameters['chord'] / 2 + 1.25 * aircraftParameters['rotorRadius'],
                          fuselageLength / 2,
                          cockpitLength + aircraftParameters['chord'] / 2,
                          cockpitLength / 2
                          ])

    return np.sum(masses*locations)/np.sum(masses)


def max_wing_loads():
    F_ay = const['gravityMars'] * (aircraftParameters['rotorMass'] + aircraftParameters['wingMass'] - const['takeOffLoad']*aircraftParameters['totalMass']/2)
    F_az = -const['gravityMars']/aircraftParameters['chord'] * (aircraftParameters['wingspan']/4 * (aircraftParameters['wingMass']-const['takeOffLoad']*aircraftParameters['totalMass']) + aircraftParameters['rotorMass'])

    F_bz = -F_az

    M_az = const['gravityMars'] * aircraftParameters['totalMass'] * aircraftParameters['cd'] / aircraftParameters['cl'] * aircraftParameters['wingspan'] / 4

    Fa = np.array([0, F_ay, F_az])
    Fb = np.array([0, F_bz, 0])
    M = np.array([0, 0, M_az])

    return Fa, Fb, M


def max_rotor_loads():
    c = aircraftParameters['rotorRadius'] / 25
    s_rotor = c * aircraftParameters['rotorRadius']
    tipSpeed = -0.88 * const['cruiseSpeed'] + 268.87
    blade_drag = aircraftParameters['cd'] * 1/2 * const['airDensity'] * tipSpeed**2 * s_rotor

    F_x = blade_drag
    F_z = (aircraftParameters['rotorMass'] - const['takeOffLoad'] * const['gravityMars'] * aircraftParameters['totalMass'] ) / 24

    M_x = aircraftParameters['rotorRadius']/2 * F_z
    M_z = -F_x * aircraftParameters['rotorRadius']/2

    return np.array([F_x, 0, F_z]), np.array([M_x, 0, M_z])


def tail_pole_loads(tailLift):
    tailDrag = tailLift * aircraftParameters['cd']/aircraftParameters['cl']
    lw = 1.25 * aircraftParameters['rotorRadius']
    lf = 1.78 * 2.5


    F_x = tailDrag
    F_z = tailLift - const['gravityMars'] * aircraftParameters['tailMass']

    M_y = (lw/2 - lf/4) * (const['gravityMars'] * aircraftParameters['tailMass'] - tailLift)

    return np.array([F_x, 0, F_z]), np.array([0, M_y, 0])
