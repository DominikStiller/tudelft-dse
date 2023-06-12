from StructureClasses import Beam, Force, xflr_forces
from rotor_sizing import y_transformation
from material_properties import materials
import vibration_toolbox as vtb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from colorama import Fore
from tqdm import tqdm


def size_rotor_blades():
    global rootBladeChord, tipBladeChord


    print(Fore.WHITE + '\n### Rotor blade sizing started ###\n')
    # Define the blade
    frontBladeChord = np.flip(np.array([26,         13,          8.66666667,  6.5,         5.2,         4.33333333,
  3.71428571,  3.25      ,  2.88888889,  2.6       ,  2.36363636,  2.16666667,
  2.        ,  1.85714286,  1.73333333,  1.625     ,  1.52941176,  1.44444444,
  1.36842105,  1.3       ,  1.23809524,  1.18181818,  1.13043478,  1.08333333,
  1.04      ,  1.        ,  0.96296296,  0.92857143,  0.89655172,  0.86666667,
  0.83870968,  0.8125    ,  0.78787879,  0.76470588,  0.74285714,  0.72222222,
  0.7027027 ,  0.68421053,  0.66666667,  0.65      ,  0.63414634,  0.61904762,
  0.60465116,  0.59090909,  0.57777778,  0.56521739,  0.55319149,  0.54166667,
  0.53061224,  0.52      ,  0.50980392,  0.5       ,  0.49056604,  0.48148148,
  0.47272727,  0.46428571,  0.45614035,  0.44827586,  0.44067797,  0.43333333,
  0.42622951,  0.41935484,  0.41269841,  0.40625   ,  0.4       ,  0.39393939,
  0.3880597 ,  0.38235294,  0.37681159,  0.37142857,  0.36619718,  0.36111111,
  0.35616438,  0.35135135,  0.34666667,  0.34210526,  0.33766234,  0.33333333,
  0.32911392,  0.325     ,  0.32098765,  0.31707317,  0.31325301,  0.30952381,
  0.30588235,  0.30232558,  0.29885057,  0.29545455,  0.29213483,  0.28888889,
  0.28571429,  0.2826087 ,  0.27956989,  0.27659574,  0.27368421,  0.27083333,
  0.26804124,  0.26530612,  0.26262626,  0.26]))
    frontBladeTwist = np.flip(np.radians(np.array([42.43378774, 41.14750545, 39.24159364, 36.96781713, 34.55063734, 32.15582014,
 29.8877177 , 27.80064951, 25.91407413, 24.22611721, 22.72353279, 21.38815179,
 20.20064245, 19.14248423, 18.19686297, 17.34896345, 16.58595097, 15.89680929,
 15.27212576, 14.70387011, 14.18518831, 13.71022028, 13.27394304, 12.87203799,
 12.50077967, 12.156943  , 11.83772601, 11.54068562, 11.26368411, 11.00484444,
 10.76251294, 10.53522794, 10.32169343, 10.12075678,  9.93138994,  9.75267339,
  9.58378249,  9.42397582,  9.27258518,  9.12900693,  8.99269456,  8.86315225,
  8.73992928,  8.62261514,  8.51083527,  8.40424736,  8.30253805,  8.20542009,
  8.11262975,  8.02392466,  7.93908175,  7.85789554,  7.78017659,  7.70575005,
  7.63445451,  7.56614085,  7.50067124,  7.43791833,  7.3777644 ,  7.32010074,
  7.26482696,  7.21185048,  7.16108605,  7.11245529,  7.06588638,  7.02131364,
  6.97867738,  6.93792356,  6.89900368,  6.86187461,  6.82649849,  6.79284272,
  6.7608799 ,  6.73058793,  6.70195011,  6.67495526,  6.64959799,  6.62587902,
  6.60380552,  6.58339167,  6.56465925,  6.5476384 ,  6.53236861,  6.51889984,
  6.50729399,  6.49762663,  6.48998925,  6.48449195,  6.48126692,  6.48047276,
  6.48230009,  6.48697869,  6.49478693,  6.50606422,  6.52122773,  6.54079533,
  6.5654176 ,  6.59592351,  6.63338738,  6.67922986])))

    rearBladeChord = np.flip(np.array([26.,         13.,          8.66666667,  6.5,         5.2,         4.33333333,
  3.71428571,  3.25      ,  2.88888889,  2.6       ,  2.36363636,  2.16666667,
  2.        ,  1.85714286,  1.73333333,  1.625     ,  1.52941176,  1.44444444,
  1.36842105,  1.3       ,  1.23809524,  1.18181818,  1.13043478,  1.08333333,
  1.04      ,  1.        ,  0.96296296,  0.92857143,  0.89655172,  0.86666667,
  0.83870968,  0.8125    ,  0.78787879,  0.76470588,  0.74285714,  0.72222222,
  0.7027027 ,  0.68421053,  0.66666667,  0.65      ,  0.63414634,  0.61904762,
  0.60465116,  0.59090909,  0.57777778,  0.56521739,  0.55319149,  0.54166667,
  0.53061224,  0.52      ,  0.50980392,  0.5       ,  0.49056604,  0.48148148,
  0.47272727,  0.46428571,  0.45614035,  0.44827586,  0.44067797,  0.43333333,
  0.42622951,  0.41935484,  0.41269841,  0.40625   ,  0.4       ,  0.39393939,
  0.3880597 ,  0.38235294,  0.37681159,  0.37142857,  0.36619718,  0.36111111,
  0.35616438,  0.35135135,  0.34666667,  0.34210526,  0.33766234,  0.33333333,
  0.32911392,  0.325     ,  0.32098765,  0.31707317,  0.31325301,  0.30952381,
  0.30588235,  0.30232558,  0.29885057,  0.29545455,  0.29213483,  0.28888889,
  0.28571429,  0.2826087 ,  0.27956989,  0.27659574,  0.27368421,  0.27083333,
  0.26804124,  0.26530612,  0.26262626,  0.26      ]))
    rearBladeTwist = np.flip(np.radians(np.array([ 31.8609133,   33.78500837,  35.43535859, -23.15916849, -25.97132754,
 -31.18604613, -40.90954726, -48.97296653, -60.27684071, -14.97638318,
 -13.60117457, -12.50470998, -11.61168686,  37.96510886,  37.37429069,
  36.74423084,  36.07711508,  35.38519305,  34.67723838,  33.96022864,
  33.23987975,  32.52087892,  31.80704483,  31.10145937,  30.40658074,
  29.72434099,  29.05622968,  28.40336557,  27.76655785,  27.14635821,
  26.5431053 ,  25.95696239,  25.38794934,  24.83596973,  24.30083364,
  23.78227696,  23.27997755,  22.7935688 ,  22.32265093,  21.8668004,
  21.4255777 ,  20.99853373,  20.58521509,  20.18516831,  19.79794333,
  19.4230962 ,  19.0601913 ,  15.17840633,  13.82929334,  12.89977478,
  12.18362701,  11.60043421,  11.10896107,  10.68483948,  10.31233707,
   9.98065934,   9.68206194,   9.41079541,   9.1624739 ,   8.93367707,
   8.72168811,   8.52431527,   8.33976668,   8.16656048,   8.00345875,
   7.8494182 ,   7.70355272,   7.56510447,   7.4334215 ,   7.30793995,
   7.18817005,   7.07368484,   6.9641111 ,   6.85912203,   6.75843132,
   6.6617884 ,   6.56897466,   6.47980041,   6.39410272,   6.31174383,
   6.23261033,   6.15661297,   6.08368736,   6.01379553,   5.94692877,
   5.88311202,   5.82241055,   5.76493984,   5.71088033,   5.66049991,
   5.61418872,   5.57251537,   5.53632148,   5.50689046,   5.48627199,
   5.47797263,   5.48865277,   5.53334028,   5.65971286,   5.64894463])))

    cutout = 0.47
    discretization = 5

    frontTwist = np.zeros(discretization * (np.size(frontBladeTwist)-1))
    for i in range(np.size(frontBladeTwist)-1):
        frontTwist[i*discretization:i*discretization+discretization] = np.linspace(frontBladeTwist[i], frontBladeTwist[i+1], discretization)

    rearTwist = np.zeros(discretization * (np.size(rearBladeTwist)-1))
    for i in range(np.size(frontBladeTwist)-1):
        rearTwist[i*discretization:i*discretization+discretization] = np.linspace(rearBladeTwist[i], rearBladeTwist[i+1], discretization)

    frontTwist = frontTwist[:np.size(frontTwist)-round(cutout*np.size(frontTwist))]
    l_front = np.linspace(-R, -R * cutout, np.size(frontTwist))

    rearTwist = rearTwist[:np.size(rearTwist)-round(cutout*np.size(rearTwist))]
    l_rear = np.linspace(-R, -R * cutout, np.size(frontTwist))

    frontChord = np.linspace(frontBladeChord[0],
                             frontBladeChord[np.size(frontBladeChord) - round(cutout * np.size(frontBladeChord))],
                             np.size(frontTwist))

    rearChord = np.linspace(rearBladeChord[0],
                             rearBladeChord[np.size(rearBladeChord) - round(cutout * np.size(rearBladeChord))],
                             np.size(rearTwist))

    frontSect = np.vstack((Airfoil['x'], Airfoil['z']))
    frontSect = np.ones((np.size(l_front), 2, 1)) * frontSect * np.reshape(frontChord, (np.size(frontChord), 1, 1))
    frontSect = y_transformation(frontTwist, frontSect)

    rearSect = np.vstack((Airfoil['x'], Airfoil['z']))
    rearSect = np.ones((np.size(l_rear), 2, 1)) * rearSect * np.reshape(rearChord, (np.size(rearChord), 1, 1))
    rearSect = y_transformation(rearTwist, rearSect)

    global Xac_rotor, Zac_rotor
    Xac_rotor = np.max(Airfoil['x']) * frontChord[-1] / 4
    Zac_rotor = 0.077 * frontChord[0] / 20

    frontBlade = Beam(
        width=frontSect[:, 0].T,
        height=frontSect[:, 1].T,
        length=l_front,
        cross_section=frontSect,
        material='CFRCy',
        fixing_points=np.array([[Xac_rotor], [Zac_rotor]]) * np.ones(np.size(l_front))
    )

    rearBlade = Beam(
        width=rearSect[:, 0].T,
        height=rearSect[:, 1].T,
        length=l_rear,
        cross_section=rearSect,
        material='CFRCy',
        fixing_points=np.array([[Xac_rotor], [Zac_rotor]]) * np.ones(np.size(l_front))
    )

    # Define the applied forces
    liftOffperBlade = np.ones(np.shape(frontTwist)) * 1.1 * MTOM * g / 32 / np.size(frontTwist)
    application = np.ones(np.shape(frontTwist)) * np.array([[Xac_rotor], [-R], [Zac_rotor]])
    application[1] = l_front

    liftOffForce_front = Force(
        magnitude=liftOffperBlade * np.array(
            [
                [-1 / 54.178],
                [0],
                [1]
            ]
        ),
        point_of_application=application
    )

    application2 = np.ones(np.shape(rearTwist)) * np.array([[Xac_rotor], [-R], [Zac_rotor]])
    application2[1] = l_rear
    liftOffForce_rear = Force(
        magnitude=liftOffperBlade * np.array(
            [
                [-1 / 54.178],
                [0],
                [1]
            ]
        ),
        point_of_application=application2
    )

    diff = 100
    rotorMass = 250 / 9 / np.size(l_front)
    print('Iteration begun')
    while diff > 0.01:
        frontBlade.unload()
        rearBlade.unload()

        rotationForce_front = rotorMass * rpm**2 * l_front
        rotatingForce_front = Force(
            magnitude=np.vstack((np.zeros(np.size(rotationForce_front)), rotationForce_front, np.zeros(np.size(rotationForce_front)))),
            point_of_application=application
        )

        rotationForce_rear = rotorMass * rpm ** 2 * l_rear
        rotatingForce_rear = Force(
            magnitude=np.vstack((np.zeros(np.size(rotationForce_rear)), rotationForce_rear,
                                 np.zeros(np.size(rotationForce_rear)))),
            point_of_application=application
        )

        frontBlade.add_loading(liftOffForce_front)
        frontBlade.add_loading(rotatingForce_front)
        frontBlade.plot_internal_loading()

        rearBlade.add_loading(liftOffForce_rear)
        rearBlade.add_loading(rotatingForce_rear)
        rearBlade.plot_internal_loading()

        frontBlade.InternalStress(0, 0, 0)
        frontBlade.calculate_mass()

        rearBlade.InternalStress(0, 0, 0)
        rearBlade.calculate_mass()

        diff = np.maximum(np.abs(rotorMass - frontBlade.m - 10) / (rotorMass - 10),
                          np.abs(rotorMass - rearBlade.m - 10) / (rotorMass - 10))
        rotorMass = np.hstack((np.sum(frontBlade.masses(), 0), np.array([0]))) + 10 * np.ones(np.size(l_front)) / np.size(l_front)

    def I(w, t):
        return w * t ** 3 + 2 * w * t * (0.073 - t) / 2 + t * (0.073 - 2 * t) ** 3 / 6

    required_I = 2.300534277610573e-06
    for i in np.linspace(0.001, frontChord[-1], 1000):
        if I(i, 0.001) > required_I:
            break

    print(f'A hollow square beam of thickness 1[mm] and width {round(i*1e3, 2)} [mm] is needed before the cutout')
    rod_weight_1 = (2 * i*0.001 + 2*0.001 * (0.0840609-2*0.001)) * 0.5*R * materials["CFRCy"].rho
    print(f'This beam weights {rod_weight_1} [kg]')

    m_r = rotorMass * 24

    print(f'Each front blade weights {np.round(frontBlade.m, 2) + 10} kg, including 10 kg of reinforcements')
    print(f'Each rear blade weights {np.round(rearBlade.m, 2) + 10} kg, including 10 kg of reinforcements')
    print(f'Total rotor mass = {np.round(12 * (frontBlade.m + 10) + 12 * (rearBlade.m + 10), 2)} kg')

    wn_rotor, x_rotor, U_rotor = rotor_vibrations(frontBlade)
    wn_rotor, x_rotor, U_rotor = rotor_vibrations(rearBlade)
    m_prop = np.sum(m_r) + rod_weight_1
    print(Fore.BLUE + f'The total mass of the propulsion subsystem is {round(m_prop, 2)} [kg]')
    return frontBlade, rearBlade, m_prop


def plot_mode_response(x, U):
    fig, ax = plt.subplots()
    ax.plot(x, U * 1e3, label=[f'Mode {i + 1}' for i in range(np.shape(U)[1])])
    ax.set_xlabel('Span [m]')
    ax.set_ylabel('Displacement [mm]')

    ax.grid(True)
    ax.legend()

    fig.tight_layout()
    plt.show()


def equivalent_load(deflection, position, E, I):
    radius_of_curvature = position**2 / deflection
    curvature = 1 / radius_of_curvature

    # Assuming linearly elastic:
    M = curvature * E * I
    P = M / position
    return P


def rotor_vibrations(rotorBlade):
    # Parameters of the rod
    cutout = 0.47
    L1 = R * cutout
    E1 = materials['CFRCy'].E
    I0 = np.mean(np.sum((rotorBlade.Bi * rotorBlade.z[:-1]**2), 0))
    A0 = np.mean(np.sum(rotorBlade.Bi, 0))

    print('\nVibration data:')
    print(f'Average moment of inertia of the airfoil = {I0}')

    bc = 2  # Clamped - free
    modes = 3

    parameters = np.array([E1, I0, materials['CFRCy'].rho, A0, L1])
    w, x, U = vtb.euler_beam_modes(n=modes, bctype=bc, beamparams=parameters)

    print(f'Original natural freq: {np.min(w)}')

    reinforcement_area = 0
    z_NA = np.sum(rotorBlade.Bi*rotorBlade.z[:-1], 0)/np.sum(rotorBlade.Bi, 0)
    z_max = np.max(rotorBlade.z, 0)
    arm = np.mean(z_max - z_NA)
    dA = 0.0001
    while np.min(w) < 20:
        A0 += dA
        I0 += dA * arm**2
        reinforcement_area += dA
        if reinforcement_area >= 1:
            raise ValueError('The added area is too large')

        parameters = np.array([E1, I0, materials['CFRCy'].rho, A0, L1])
        w, x, U = vtb.euler_beam_modes(n=modes, bctype=bc, beamparams=parameters)

    print(f'Natural frequencies = {w} [Hz]')
    print(f'Maximum deflection = {np.max(np.abs(U))} [m]')
    print(f'Required reinforcement area = {reinforcement_area}')
    print(f'Required reinforcement mass = {reinforcement_area*R*(1-cutout)*materials["CFRCy"].rho} kg')
    plot_mode_response(x, U)

    # Calculate equivalent load and additional stress
    max_deflection = U[np.where(np.abs(U) == np.max(np.abs(U)))]
    max_deflection_pos = x[np.where(np.abs(U) == np.max(np.abs(U)))[0]]
    P = equivalent_load(deflection=max_deflection,
                        position=max_deflection_pos,
                        E=E1,
                        I=I0)

    print(f'Equivalent load due to vibrations in rotor blades is {P} [N]')
    return w, x, U


def wing_vibrations():
    # Define parameters
    L = span[best] / 2
    I = np.mean(np.sum((best_wing.Bi * best_wing.z[:-1] ** 2), 0))
    A = np.mean(np.sum(best_wing.Bi, 0))
    parameters = np.array([materials['CFRPeek'].E, I, materials['CFRPeek'].rho, A, L])

    # Define analysis conditions
    bc = 6  # Clamped-pinned
    modes = 3

    # Perform analysis
    w, x, U = vtb.euler_beam_modes(n=modes, bctype=bc, beamparams=parameters)

    # Report results
    print(Fore.WHITE + '\nVibration data:')
    print(f'Average moment of inertia of the airfoil = {I}')
    print(f'Natural frequency of the wing = {w}')
    print(f'Max deflection = {np.max(np.abs(U))} [m]')
    plot_mode_response(x, U)

    return w, x, U


def size_wing(span, chord_root, taper, wing_model=None):
    l = np.linspace(-span, 0, 100)
    chord_array = np.linspace(chord_root*taper, chord_root, np.size(l))

    # Define the geometry
    global Xac_wing, Zac_wing
    Xac_wing = np.max(Airfoil['x']) * chord_array[-1] / 4
    Zac_wing = 0.077 * chord_array[-1]


    section = np.vstack((Airfoil["x"], Airfoil["z"])) * np.reshape(np.vstack((chord_array, chord_array)), (np.size(chord_array), 2, 1))

    wing = Beam(
        width=np.array([Airfoil['x']]).T * chord_array,
        height=np.array([Airfoil['z']]).T * chord_array,
        length=l,
        cross_section=section,
        material='CFRPeek',
        fixing_points=np.array([[Xac_wing], [Zac_wing]]) * np.ones(np.size(l))
    )

    theta = np.arctan(fuselage_height / span)

    # Define the forces during TO
    bracing_TO_mag = g / np.sin(theta) * (1.1 * MTOM / 2 - (mr/2 + m_e))
    bracing_TO = Force(
        magnitude=bracing_TO_mag * np.array(
            [
                [0],
                [np.cos(theta)],
                [-np.sin(theta)]
            ]
        ),
        point_of_application=np.array(
            [
                [Xac_wing],
                [-span],
                [Zac_wing]
            ]
        )
    )
    engine_and_rotor_weight = Force(
        magnitude=np.array(
            [
                [0],
                [0],
                [-(mr/2 + m_e) * g]
            ]
        ),
        point_of_application=np.array(
            [
                [Xac_wing],
                [-span],
                [Zac_wing]
            ]
        )
    )
    liftOffLoad = Force(
        magnitude=np.array(
            [
                [0],
                [0],
                [1.1 * MTOM * g / 2]
            ]
        ),
        point_of_application=np.array(
            [
                [Xac_wing],
                [-span],
                [Zac_wing]
            ]
        )
    )
    cruiseThrust = Force(
        magnitude=np.array(
            [
                [790 / 2],
                [0],
                [0]
            ]
        ),
        point_of_application=np.array(
            [
                [Xac_wing],
                [-span],
                [Zac_wing]
            ]
        )
    )

    # Define loads during cruise
    if wing_model is None:
        aerodynamic_forces = xflr_forces('Test_xflr5_file.csv', q, float(span))
    else:
        cl, cd, y_L = xflr_forces('wings.csv', q, float(span), adrian=wing_model)
        dy = np.abs(wing.y[:-1] - wing.y[1:])
        c_arr = (chord_array[:-1] + chord_array[1:]) / 2

        Cl = np.zeros(np.size(dy))
        Cd = np.zeros(np.size(dy))
        app = np.zeros(np.size(dy))
        for i in range(np.size(dy)):
            indx = (np.abs(wing.y[i] - y_L)).argmin()
            Cl[i] = cl[indx]
            Cd[i] = cd[indx]
            app[i] = wing.y[i]

        lift = Cl * q * c_arr * dy
        drag = Cd * q * c_arr * dy

        aerodynamic_forces = Force(
            magnitude=np.vstack((-drag, np.zeros(np.size(lift)), lift)),
            point_of_application=np.vstack((Xac_wing * np.ones(np.size(lift)),
                                            app,
                                            Zac_wing * np.ones(np.size(lift))))
        )

    liftMoment = -np.dot(aerodynamic_forces.F[2], aerodynamic_forces.application[1])
    bracing_cr_mag = 1 / (span * np.sin(theta)) * (liftMoment - span * g * (mr/2 + m_e))
    bracing_cruise = Force(
        magnitude=bracing_cr_mag * np.array(
            [
                [0],
                [np.cos(theta)],
                [-np.sin(theta)]
            ]
        ),
        point_of_application=np.array(
            [
                [Xac_wing],
                [-span],
                [Zac_wing]
            ]
        )
    )

    # Apply loads and size
    wing.add_loading(engine_and_rotor_weight)
    wing.add_loading(cruiseThrust)
    wing.add_loading(bracing_cruise)
    wing.add_loading(aerodynamic_forces)
    wing.plot_internal_loading()
    wing.InternalStress(0, 0, 0)
    f_loading_Cr = wing.f_loading
    thickness2 = wing.t
    B2 = wing.Bi

    # Apply loads and size
    wing.unload()
    wing.add_loading(liftOffLoad)
    wing.add_loading(engine_and_rotor_weight)
    wing.add_loading(bracing_TO)
    wing.plot_internal_loading()
    wing.InternalStress(0, 0, 0)
    f_loading_TO = wing.f_loading
    thickness = wing.t
    B1 = wing.Bi

    # Choose the most critical structure
    wing.t = np.maximum(thickness, thickness2)
    wing.Bi = np.maximum(B1, B2)
    f_loading_abs = np.maximum(np.abs(f_loading_TO), np.abs(f_loading_Cr))

    wing.calculate_mass()
    print(f'Mass of each wing = {np.round(wing.m, 2)}kg')
    return wing, f_loading_abs


def size_tail():
    # span, 14.23/2  chord 2.7106, tip = 1.0834 NACA0012

    # Assumptions
    m = 25
    tail_to_wing_lift = 0.1
    lift_to_drag_tail = 15
    extra_force = 500  # For control
    tail_taper = np.linspace(1, 0.7, m)
    vTailTaper = tail_taper

    ### Horizontal stabilizer ###
    # Define the geometry
    tailSpan = 4.19
    tailChord = np.linspace(1.1996, 2.998, m)
    NACA0012 = pd.read_csv('NACA 0012.dat', delimiter="\s+", dtype=float, skiprows=1, names=["x", "z"])

    l = np.linspace(-tailSpan, 0, m)
    x = np.reshape(NACA0012['x'], (len(NACA0012['x']), 1)) * tailChord
    z = np.reshape(NACA0012['z'], (len(NACA0012['z']), 1)) * tailChord
    section = np.vstack((NACA0012['x'], NACA0012['z'])) * np.reshape(np.vstack((tailChord, tailChord)), (m, 2, 1))

    Xac = np.max(NACA0012['x']) * tailChord[0] / 4
    Zac = 0.077 * tailChord[0]

    hStabilizer = Beam(
        width=x,
        height=z,
        length=l,
        cross_section=section,
        material='Al/Si',
        fixing_points=np.array([[Xac], [Zac]]) * np.ones(m)
    )

    # Define the loads
    lift = Force(  # Assuming uniform distribution for now, will change later
        magnitude=(tail_to_wing_lift*MTOM*g + extra_force) / m * np.vstack((np.zeros((2, m)), np.ones(m))),
        point_of_application=np.vstack((Xac*np.ones(m), l, Zac*np.ones(m)))
    )

    drag = Force(  # Assuming uniform distribution for now, will change later
        magnitude=(tail_to_wing_lift*MTOM*g + extra_force) / (m * lift_to_drag_tail) * np.vstack((-1*np.ones(m), np.zeros((2, m)))),
        point_of_application=np.vstack((Xac*np.ones(m), l, Zac*np.ones(m)))
    )

    # Apply loads and size
    hStabilizer.add_loading(lift)
    hStabilizer.add_loading(drag)
    hStabilizer.plot_internal_loading()
    hStabilizer.InternalStress(0, 0, 0)
    hStabilizer.calculate_mass()
    print(f"Horizontal stabilizer's mass is {2*np.round(hStabilizer.m)} [kg]")

    ### Vertical stabilizer ###
    # Define geometry

    vTailSpan = 14.23/2
    vTailChord = np.linspace(1.0834, 2.7106, m)


    lv = np.linspace(-vTailSpan, 0, m)
    xv = np.reshape(NACA0012['x'], (len(NACA0012['x']), 1)) * vTailChord
    zv = np.reshape(NACA0012['z'], (len(NACA0012['x']), 1)) * vTailChord
    section = np.vstack((NACA0012['x'], NACA0012['z'])) * np.reshape(np.vstack((vTailChord, vTailChord)), (m, 2, 1))

    Xac = np.max(NACA0012['x']) * vTailChord[0] / 4
    Zac = 0.077 * vTailChord[0]

    vStabilizer = Beam(
        width=xv * vTailTaper,
        height=zv * vTailTaper,
        length=lv,
        cross_section=section,
        material='Al/Si',
        fixing_points=np.array([[Xac], [Zac]]) * np.ones(m)
    )

    # Define loads
    restoring = Force(
        magnitude=extra_force / m * np.vstack((-np.ones(m)/lift_to_drag_tail, np.zeros(m), np.ones(m))),
        point_of_application=np.vstack(
            ((Xac * np.min(vTailTaper) + (1-np.min(vTailTaper))*vTailChord) * np.ones(m),
             lv,
             Zac * np.min(vTailTaper) * np.ones(m))
        )
    )

    # Apply loads and size
    vStabilizer.add_loading(restoring)
    vStabilizer.plot_internal_loading()
    vStabilizer.InternalStress(0, 0, 0)
    vStabilizer.calculate_mass()
    print(f"Vertical stabilizer's mass is {np.round(vStabilizer.m)} [kg]")


    ### Tail pole ###
    margin = 50  # [kg] Mass margin for the tail assembly and reinforcements
    boomMoments = vStabilizer.m_loading[-1] + lengthTailPole * (np.array(
        [
            [-hStabilizer.f_loading[-1][2][0]],
            [0],
            [0]
        ]
    ) + np.array(
        [
            [0],
            [0],
            [-vStabilizer.f_loading[-1][2][0]]
        ]
    ) + np.array(
        [
            [g * (hStabilizer.m + vStabilizer.m + margin)],
            [0],
            [0]
        ]
    ))

    s_max = materials['CFRP'].compressive
    t_min = 0.001
    Mx = boomMoments[0]
    Mz = boomMoments[2]
    Fy = hStabilizer.f_loading[-1][0][0] + vStabilizer.f_loading[-1][1][0]
    theta = np.arctan(Mx/Mz)
    D = (Fy + np.sqrt(Fy**2 + 8 * s_max * np.pi * t_min * np.abs(Mz*np.cos(theta) + Mx*np.sin(theta)))) / (2 * s_max * np.pi * t_min)

    tailPoleMass = np.pi * D * 0.001 * R * materials['CFRPeek'].rho

    print(f'Tail pole mass = {tailPoleMass} [kg]')
    print(Fore.BLUE + f'The total tail group mass is {2*hStabilizer.m + vStabilizer.m + tailPoleMass + margin} [kg], '
          f'including {margin} [kg] of margin')
    return hStabilizer, vStabilizer, tailPoleMass


def size_body(fuselage_height=1.67, cabin_length=2, full_length=6.15):
    r = fuselage_height/2
    aft_cone_length = full_length - cabin_length - rootChord  # [m], assumed
    mat = materials['CFRPeek']
    t = 0.001
    margin = 200

    # From https://dr.ntu.edu.sg/bitstream/10356/152917/2/Damage%20Severity%20Prediction%20of%20Helicopter%20Windshield.pdf
    windshieldMaterial = materials['Polycarbonate']
    windshieldThickness = 0.005

    cabin_SA = np.pi/6 * (r/cabin_length**2) * ((r**2 + 4*cabin_length**2)**1.5 - r**3)
    main_body_SA = 2 * np.pi * r * rootChord
    aft_connection_SA = np.pi/6 * (r/aft_cone_length**2) * ((r**2 + 4*aft_cone_length**2)**1.5 - r**3)

    cabinMass = cabin_SA * windshieldThickness * windshieldMaterial.rho
    bodyMass = main_body_SA * t * mat.rho
    aftConeMass = aft_connection_SA * t * mat.rho

    bodySkinMass = bodyMass + aftConeMass + cabinMass + margin

    print(f'Cabin mass = {cabin_SA * windshieldThickness * windshieldMaterial.rho} [kg]')
    print(f'Body mass = {(main_body_SA + aft_connection_SA) * t * mat.rho} [kg]')
    print(f'Margin = {margin} [kg]')

    # Moments of inertia
    ycg = (cabinMass * (full_length - 2/3*cabin_length) + bodyMass * (full_length - cabin_length - rootChord/2) +
           aftConeMass * 2/3 * (full_length-cabin_length-rootChord)) / (cabinMass + bodyMass + aftConeMass)
    print(f'Body cg is {full_length - ycg} [m] behind the nose')

    # Cabin
    Ix_c = 2 / 3 * cabinMass * r * windshieldThickness
    Iy_c = Iz_c = cabinMass * (r*windshieldThickness / 3 + cabin_length*windshieldThickness/9)

    # Main body
    Ix_b = bodyMass * r * t
    Iy_b = Iz_b = bodyMass * r * t / 2

    # Aft cone
    Ix_a = 2 / 3 * aftConeMass * r * t
    Iy_a = Iz_a = aftConeMass * (r*t / 3 + (full_length-cabin_length-rootChord)*t/9)

    # Combine
    Ix = Ix_c + cabinMass * ((full_length - 2/3*cabin_length) - ycg)**2 + \
         Ix_b + bodyMass * ((full_length - cabin_length - rootChord/2) - ycg)**2 + \
         Ix_a + aftConeMass * (2/3 * (full_length-cabin_length-rootChord) - ycg)**2

    Iy = Iy_c + Iy_b + Iy_a
    Iz = Iz_c + Iz_b + Iz_a

    return bodySkinMass, Ix, Iy, Iz



if __name__ == '__main__':
    global Xac_rotor, Zac_rotor, Xac_wing, Zac_wing, rootBladeChord, tipBladeChord

    # Constants
    R = 10.4
    MTOM = 3000
    g = 3.71
    rpm = 183 * 2 * np.pi / 60
    q = 0.5 * 0.01 * 112 ** 2
    rootChord = 3.35
    b = 16.8
    fuselage_height = 1.67
    m_e = 50

    cabin_length = 1.73
    x_ac_w = cabin_length + rootChord / 4
    x_cg = x_ac_w + 1
    x_ac_t = x_cg + 10
    lengthTailPole = x_ac_t - cabin_length - rootChord

    Airfoil = pd.read_csv("S1223.dat", delimiter="\s+", dtype=float, skiprows=1, names=["x", "z"])

    # Rotors
    frontBlade, rearBlade, mr = size_rotor_blades()

    # Wings
    print(Fore.WHITE + '\n### Wing sizing started ###\n')
    span = np.flip(np.array([45, 40, 35, 30]))
    rootChord = np.flip(np.array([4, 4.3333, 4.9481, 5.78]))
    tr = 0.5
    min_mass = 1e6
    for i in range(len(span)):
        b, c = span[i], rootChord[i]
        wing, F = size_wing(span=b/2, chord_root=c, taper=tr, wing_model=i)
        wing.calculate_mass()

        if wing.m < min_mass:
            best = i
            best_wing = wing
            min_mass = wing.m
            best_wing.f_loading = F

        print(f'Configuration {i} explored, m = {wing.m} [kg]')

    print(Fore.BLUE + f'Configuration {best} is the best, with a total mass of {2*round(min_mass, 2)} [kg] per wing')
    wn_wing, x_wing, U_wing = wing_vibrations()
    n_rivets_0 = best_wing.design_joint(b=span[best]/2)
    n_rivets_1 = best_wing.design_joint(b=0)

    print(Fore.WHITE + '\n### Tail sizing started ###\n')
    hStabilizer, vStabilizer, tailPoleMass = size_tail()
