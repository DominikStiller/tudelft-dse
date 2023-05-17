from scipy.interpolate import InterpolatedUnivariateSpline
from .constants import const, aircraftParameters
import numpy as np


def airfoilArea():
    R = aircraftParameters['rotorRadius']
    x_cord_top = np.flip(
        [1, 0.99838, 0.99417, 0.98825, 0.98075, 0.97111, 0.95884, 0.94389, 0.92639, 0.90641, 0.88406, 0.85947, 0.83277,
         0.80412, 0.77369, 0.74166, 0.70823, 0.6736, 0.63798, 0.60158, 0.56465, 0.52744, 0.49025, 0.4534, 0.41721,
         0.38193,
         0.34777, 0.31488, 0.28347, 0.2537, 0.22541, 0.19846, 0.17286, 0.14863, 0.12591, 0.10482, 0.08545, 0.06789,
         0.05223,
         0.03855, 0.02694, 0.01755, 0.01028, 0.00495, 0.00155, 0.00005])
    x_chord_bottom = np.array([0.00005, 0.00044, 0.00264, 0.00789, 0.01718, 0.03006, 0.04627,
                      0.06561, 0.08787, 0.11282, 0.1402, 0.17006, 0.20278, 0.2384, 0.27673, 0.3175, 0.36044, 0.40519,
                      0.45139, 0.4986, 0.54639,
                      0.59428, 0.64176, 0.68832, 0.73344, 0.7766, 0.81729, 0.855, 0.88928, 0.91966, 0.94573, 0.96693,
                      0.98255, 0.99268, 0.99825, 1])
    y_cord_top = np.flip(
        [0, 0.00126, 0.00494, 0.01037, 0.01646, 0.0225, 0.02853, 0.03476, 0.04116, 0.04768, 0.05427, 0.06089, 0.06749,
         0.07402, 0.08044, 0.08671, 0.09277, 0.09859, 0.10412, 0.10935, 0.11425, 0.11881, 0.12303, 0.12683, 0.13011,
         0.13271, 0.13447, 0.13526, 0.13505, 0.13346, 0.13037, 0.12594, 0.12026, 0.11355, 0.10598, 0.0977, 0.08879,
         0.0794, 0.06965, 0.05968, 0.04966, 0.03961, 0.02954, 0.01969, 0.01033, 0.00178])
    y_chord_bottom = np.array([0.00178, -0.00561, -0.0112, -0.01427, -0.0155, -0.01584, -0.01532, -0.01404, -0.01202, -0.00925,
                      -0.00563, -0.00075, 0.00535, 0.01213, 0.01928, 0.02652, 0.03358, 0.04021, 0.04618, 0.05129,
                      0.05534, 0.0582, 0.05976, 0.05994, 0.05872, 0.05612, 0.05219, 0.04706, 0.04088, 0.03387, 0.02624,
                      0.01822, 0.0106, 0.00468, 0.00115, 0])

    S1223_top = InterpolatedUnivariateSpline(x_cord_top, y_cord_top)
    S1223_bottom = InterpolatedUnivariateSpline(x_chord_bottom, y_chord_bottom)

    Area_top = S1223_top.integral(0, 1)
    Area_bot = S1223_bottom.integral(0, 1)
    Area = (Area_top - Area_bot) * R / 20

    length = np.sum(np.sqrt((x_cord_top[1:]-x_cord_top[:-1])**2 + (y_cord_top[1:]-y_cord_top[:-1])**2)) + \
             np.sum(np.sqrt((x_chord_bottom[1:] - x_chord_bottom[:-1]) ** 2 + (y_chord_bottom[1:] - y_chord_bottom[:-1]) ** 2))

    return Area, length


def calculate_cg():
    cockpitLength = 1.78
    fuselageLength = cockpitLength + aircraftParameters['chord']

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


def max_rotor_loads(airfoil_shape):
    initial_fillfactor = 0.1
    initial_area, length = airfoilArea()
    initial_mass = const['bladeDensity'] * aircraftParameters['rotorRadius'] * initial_area * initial_fillfactor

    c = aircraftParameters['rotorRadius'] / 20
    s_rotor = c * aircraftParameters['rotorRadius']
    tipSpeed = -0.88 * const['cruiseSpeed'] + 268.87
    blade_drag = const['cd'] * 1/2 * const['airDensity'] * tipSpeed**2 * s_rotor
    omega = tipSpeed / aircraftParameters['rotorRadius']

    F_x = F_x_init = blade_drag
    F_y = F_y_init = (initial_mass) * aircraftParameters['rotorRadius'] * omega**2 / 2
    F_z = F_z_init = initial_mass*const['gravityMars'] -( const['takeOffLoad'] * const['gravityMars'] * const['maxMass'] ) / 24

    Mx = Mx_init = aircraftParameters['rotorRadius']/2 * F_z_init
    Mz = Mz_init = -F_x_init * aircraftParameters['rotorRadius']/2

    x_av, z_av = np.mean(airfoil_shape[0]), np.mean(airfoil_shape[1])

    ite = 0
    fill_factor = 0
    max_stress = 10e1000
    while np.max(np.abs(max_stress)) > const['allowedStress']:
        fill_factor += 0.0001
        guess_Rotor_mass = const['bladeDensity'] * aircraftParameters['rotorRadius'] * initial_area * fill_factor
        area = fill_factor*initial_area
        ite += 1
        F_x = blade_drag
        F_y = (guess_Rotor_mass) * aircraftParameters['rotorRadius'] * omega ** 2 / 2
        F_z = guess_Rotor_mass*const['gravityMars'] -( const['takeOffLoad'] * const['gravityMars'] * const['maxMass']/0.88) / 24

        Mx = aircraftParameters['rotorRadius'] / 2 * F_z
        Mz = -F_x * aircraftParameters['rotorRadius'] / 2
        # TODO: Fix these such that they produce the AREA moment of inertia [m^4]
        Ixx = (area / np.size(airfoil_shape[0])) * np.sum((airfoil_shape[0] - x_av) ** 2)
        Izz = (area / np.size(airfoil_shape[0])) * np.sum((airfoil_shape[1] - z_av) ** 2)
        Ixz = (area / np.size(airfoil_shape[0])) * np.sum((airfoil_shape[0] - x_av) * (airfoil_shape[1] - z_av))

        max_stress = ((Mx * Izz - Mz * Ixz) * (airfoil_shape[1]-z_av) + (Mz * Ixx - Mx * Ixz) * (airfoil_shape[0])-x_av) / (
                    Ixx * Izz - Ixz ** 2) + F_y/area

    print(f'Fill factor: {fill_factor}')
    print(f'Final Rotor Mass possible: {const["bladeDensity"] * aircraftParameters["rotorRadius"] * initial_area * (fill_factor)}')
    print(f'Thickness: {(initial_area*fill_factor)/length}')
    print(f'Maximum absolute stress: {np.max(np.abs(max_stress))}')
    return np.array([F_x, F_y, F_z]), np.array([Mx, 0, Mz]), max_stress, fill_factor


def tail_pole_loads(tailLift):
    tailDrag = tailLift * aircraftParameters['cd']/aircraftParameters['cl']
    lw = 1.25 * aircraftParameters['rotorRadius']
    lf = 1.78 * 2.5


    F_x = tailDrag
    F_z = tailLift - const['gravityMars'] * aircraftParameters['tailMass']

    M_y = (lw/2 - lf/4) * (const['gravityMars'] * aircraftParameters['tailMass'] - tailLift)

    return np.array([F_x, 0, F_z]), np.array([0, M_y, 0])



def assembly_volume():
    cockpitLength = 1.78
    fuselageWidth = 1.5

    fuselageVolume = 1/3*cockpitLength*np.pi*(fuselageWidth/2)**2 + 2*np.pi*fuselageWidth/2 * aircraftParameters['chord']

    x_cord_top = np.flip(
        [1, 0.99838, 0.99417, 0.98825, 0.98075, 0.97111, 0.95884, 0.94389, 0.92639, 0.90641, 0.88406, 0.85947, 0.83277,
         0.80412, 0.77369, 0.74166, 0.70823, 0.6736, 0.63798, 0.60158, 0.56465, 0.52744, 0.49025, 0.4534, 0.41721,
         0.38193,
         0.34777, 0.31488, 0.28347, 0.2537, 0.22541, 0.19846, 0.17286, 0.14863, 0.12591, 0.10482, 0.08545, 0.06789,
         0.05223,
         0.03855, 0.02694, 0.01755, 0.01028, 0.00495, 0.00155, 0.00005])
    x_chord_bottom = [0.00005, 0.00044, 0.00264, 0.00789, 0.01718, 0.03006, 0.04627,
                      0.06561, 0.08787, 0.11282, 0.1402, 0.17006, 0.20278, 0.2384, 0.27673, 0.3175, 0.36044, 0.40519,
                      0.45139, 0.4986, 0.54639,
                      0.59428, 0.64176, 0.68832, 0.73344, 0.7766, 0.81729, 0.855, 0.88928, 0.91966, 0.94573, 0.96693,
                      0.98255, 0.99268, 0.99825, 1]
    y_cord_top = np.flip(
        [0, 0.00126, 0.00494, 0.01037, 0.01646, 0.0225, 0.02853, 0.03476, 0.04116, 0.04768, 0.05427, 0.06089, 0.06749,
         0.07402, 0.08044, 0.08671, 0.09277, 0.09859, 0.10412, 0.10935, 0.11425, 0.11881, 0.12303, 0.12683, 0.13011,
         0.13271, 0.13447, 0.13526, 0.13505, 0.13346, 0.13037, 0.12594, 0.12026, 0.11355, 0.10598, 0.0977, 0.08879,
         0.0794, 0.06965, 0.05968, 0.04966, 0.03961, 0.02954, 0.01969, 0.01033, 0.00178])
    y_chord_bottom = [0.00178, -0.00561, -0.0112, -0.01427, -0.0155, -0.01584, -0.01532, -0.01404, -0.01202, -0.00925,
                      -0.00563, -0.00075, 0.00535, 0.01213, 0.01928, 0.02652, 0.03358, 0.04021, 0.04618, 0.05129,
                      0.05534, 0.0582, 0.05976, 0.05994, 0.05872, 0.05612, 0.05219, 0.04706, 0.04088, 0.03387, 0.02624,
                      0.01822, 0.0106, 0.00468, 0.00115, 0]

    S1223_top = InterpolatedUnivariateSpline(x_cord_top, y_cord_top)
    S1223_bottom = InterpolatedUnivariateSpline(x_chord_bottom, y_chord_bottom)

    Area_top = S1223_top.integral(0, 1)
    Area_bot = S1223_bottom.integral(0, 1)
    Area = (Area_top - Area_bot) * aircraftParameters['chord']
    wingVolume = Area * aircraftParameters['wingspan']

    rotorVolume = aircraftParameters['totalRotors'] * (Area_top - Area_bot) * (aircraftParameters['rotorRadius'])**2/20

    # Assuming tail and wing have the same density
    tailVolume = aircraftParameters['tailMass'] / aircraftParameters['wingMass'] * wingVolume

    return fuselageVolume + wingVolume + rotorVolume + tailVolume
