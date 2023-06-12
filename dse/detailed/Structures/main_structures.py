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
    frontBladeChord = np.flip(np.array([0.75,       0.75,       0.75,       0.75,       0.75,       0.75,
 0.75      , 0.75      , 0.75      , 0.75      , 0.75      , 0.75,
 0.75      , 0.75      , 0.75      , 0.75      , 0.75      , 0.75,
 0.75      , 0.75      , 0.75      , 0.75      , 0.75      , 0.75,
 0.75      , 0.75      , 0.75      , 0.75      , 0.75      , 0.75,
 0.75      , 0.75      , 0.75      , 0.75      , 0.75      , 0.75,
 0.73873874, 0.71929825, 0.7008547 , 0.68333333, 0.66666667, 0.65079365,
 0.63565891, 0.62121212, 0.60740741, 0.5942029 , 0.58156028, 0.56944444,
 0.55782313, 0.54666667, 0.53594771, 0.52564103, 0.51572327, 0.50617284,
 0.4969697 , 0.48809524, 0.47953216, 0.47126437, 0.46327684, 0.45555556,
 0.44808743, 0.44086022, 0.43386243, 0.42708333, 0.42051282, 0.41414141,
 0.4079602 , 0.40196078, 0.39613527, 0.39047619, 0.38497653, 0.37962963,
 0.37442922, 0.36936937, 0.36444444, 0.35964912, 0.35497835, 0.35042735,
 0.34599156, 0.34166667, 0.33744856, 0.33333333, 0.32931727, 0.32539683,
 0.32156863, 0.31782946, 0.31417625, 0.31060606, 0.3071161 , 0.3037037,
 0.3003663 , 0.29710145, 0.29390681, 0.29078014, 0.2877193 , 0.28472222,
 0.28178694, 0.27891156, 0.27609428, 0.27333333]))
    frontBladeTwist = np.flip(np.radians(np.array([33.99075624, 29.02635514, 25.79106317, 23.48785964, 21.74919469, 20.38071012,
 19.26947833, 18.34513154, 17.56138065, 16.88643083, 16.29765888, 15.77848867,
 15.31647064, 14.9020547 , 14.52777995, 14.18772396, 13.87711921, 13.59207983,
 13.32940333, 13.08642398, 12.8609029 , 12.65094437, 12.45493128, 12.27147489,
 12.09937516, 11.93758919, 11.78520591, 11.64142554, 11.50554285, 11.37693341,
 11.25504221, 11.13937418, 11.02948621, 10.92498055, 10.82549909, 10.73071857,
 10.58391701, 10.40084922, 10.22734367, 10.06272334,  9.9063742 ,  9.75773828,
  9.6163076 ,  9.48161888,  9.35324889,  9.23081039,  9.11394851,  9.00233764,
  8.89567859,  8.79369613,  8.69613679,  8.60276691,  8.51337091,  8.42774973,
  8.34571948,  8.26711019,  8.1917647 ,  8.11953772,  8.05029491,  7.98391214,
  7.92027477,  7.85927704,  7.80082153,  7.74481867,  7.69118633,  7.63984947,
  7.5907398 ,  7.54379556,  7.49896128,  7.45618765,  7.41543143,  7.37665536,
  7.33982817,  7.30492467,  7.27192585,  7.24081906,  7.21159832,  7.18426462,
  7.15882643,  7.13530023,  7.11371129,  7.0940945 ,  7.07649547,  7.06097192,
  7.04759529,  7.03645278,  7.0276499 ,  7.02131364,  7.01759641,  7.01668105,
  7.01878726,  7.0241799 ,  7.03317974,  7.04617788,  7.06365497,  7.08620763,
  7.11458529,  7.14974278,  7.19291733,  7.24574486])))

    rearBladeChord = np.flip(np.array([0.75,       0.75,       0.75,       0.75,       0.75,       0.75,
 0.75      , 0.75      , 0.75      , 0.75      , 0.75      , 0.75,
 0.75      , 0.75      , 0.75      , 0.75      , 0.75      , 0.75,
 0.75      , 0.75      , 0.75      , 0.75      , 0.75      , 0.75,
 0.75      , 0.75      , 0.75      , 0.75      , 0.75      , 0.75,
 0.75      , 0.75      , 0.75      , 0.75      , 0.75      , 0.75,
 0.73873874, 0.71929825, 0.7008547 , 0.68333333, 0.66666667, 0.65079365,
 0.63565891, 0.62121212, 0.60740741, 0.5942029 , 0.58156028, 0.56944444,
 0.55782313, 0.54666667, 0.53594771, 0.52564103, 0.51572327, 0.50617284,
 0.4969697 , 0.48809524, 0.47953216, 0.47126437, 0.46327684, 0.45555556,
 0.44808743, 0.44086022, 0.43386243, 0.42708333, 0.42051282, 0.41414141,
 0.4079602 , 0.40196078, 0.39613527, 0.39047619, 0.38497653, 0.37962963,
 0.37442922, 0.36936937, 0.36444444, 0.35964912, 0.35497835, 0.35042735,
 0.34599156, 0.34166667, 0.33744856, 0.33333333, 0.32931727, 0.32539683,
 0.32156863, 0.31782946, 0.31417625, 0.31060606, 0.3071161 , 0.3037037,
 0.3003663 , 0.29710145, 0.29390681, 0.29078014, 0.2877193 , 0.28472222,
 0.28178694, 0.27891156, 0.27609428, 0.27333333]))
    rearBladeTwist = np.flip(np.radians(np.array([-6.42540236, 66.75500062, 63.4744395,  60.65207818, 58.14595431, 55.87845537,
 53.80144934, 51.88248394, 50.09831685, 48.43151658, 46.86852311, 45.39846971,
 44.01243121, 42.7029248 , 41.46356743, 30.46328578, 26.52970845, 23.88564539,
 21.91245802, 20.35795523, 19.08961126, 18.02826272, 17.12272754, 16.33805522,
 15.64939031, 15.03848704, 14.49160329, 13.99816545, 13.54988884, 13.14017971,
 12.76371771, 12.41615764, 12.09391215, 11.79399069, 11.51387821, 11.25144265,
 11.09342042, 11.00741149, 10.91889559, 10.82821403, 10.73569621, 10.64165453,
 10.54638152, 10.45014846, 10.35320494, 10.25577921, 10.15807895, 10.06029224,
  9.96258884,  9.86512146,  9.76802703,  9.67142804,  9.57543374,  9.4801414,
  9.38563739,  9.29199823,  9.19929164,  9.10757736,  9.01690811,  8.92733025,
  8.83888464,  8.7516072 ,  8.6655296 ,  8.58067984,  8.49708278,  8.41476068,
  8.33373367,  8.25402024,  8.17563774,  8.09860278,  8.02293175,  7.9486413,
  7.87574882,  7.80427304,  7.73423457,  7.66565665,  7.5985659 ,  7.53299321,
  7.46897489,  7.40655396,  7.34578181,  7.28672024,  7.22944414,  7.17404485,
  7.12063473,  7.06935319,  7.02037506,  6.97392224,  6.93028057,  6.88982467,
  6.85305605,  6.8206637 ,  6.79362538,  6.77338704,  6.76220626,  6.76387969,
  6.78552292,  6.84301442,  6.9862122 ,  6.97621516])))

    cutout = 0.15
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
    frontRotorMass = rotorMass
    rearRotorMass = rotorMass
    while diff > 0.01:
        frontBlade.unload()
        rearBlade.unload()

        rotationForce_front = frontRotorMass * rpm**2 * l_front
        rotatingForce_front = Force(
            magnitude=np.vstack((np.zeros(np.size(rotationForce_front)), rotationForce_front, np.zeros(np.size(rotationForce_front)))),
            point_of_application=application
        )

        rotationForce_rear = rearRotorMass * rpm ** 2 * l_rear
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

        diff = np.maximum(np.abs(rotorMass - frontBlade.m - 5) / (rotorMass - 5),
                          np.abs(rotorMass - rearBlade.m - 5) / (rotorMass - 5))
        frontRotorMass = np.hstack((np.sum(frontBlade.masses(), 0), np.array([0]))) + 5 * np.ones(np.size(l_front)) / np.size(l_front)
        rearRotorMass = np.hstack((np.sum(rearBlade.masses(), 0), np.array([0]))) + 5 * np.ones(
            np.size(l_rear)) / np.size(l_rear)

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

    print(f'Each front blade weights {np.round(frontBlade.m, 2) + 5} kg, including 5 kg of reinforcements')
    print(f'Each rear blade weights {np.round(rearBlade.m, 2) + 5} kg, including 5 kg of reinforcements')
    print(f'Total rotor mass = {np.round(12 * (frontBlade.m + 5) + 12 * (rearBlade.m + 5), 2)} kg')

    wn_rotor, x_rotor, U_rotor = rotor_vibrations(frontBlade)
    wn_rotor, x_rotor, U_rotor = rotor_vibrations(rearBlade)
    m_prop = 12 * (frontBlade.m + 5) + 12 * (rearBlade.m + 5) + rod_weight_1
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
    L1 = R
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
    L = b / 2
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

    section = np.vstack((Airfoil["x"], Airfoil["z"])) * np.reshape(chord_array, (np.size(chord_array), 1, 1))

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

    cruiseThrust = Force(
        magnitude=np.array(
            [
                [-np.sum(aerodynamic_forces.F[0])],
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
    margin = 0  # [kg] Mass margin for the tail assembly and reinforcements
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
    R = 8.2
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

    # b = 42.69
    # rootChord = 3.7948
    # tipChord = 1.8974
    # best_wing, F = size_wing(span=b/2, chord_root=rootChord, taper=rootChord/tipChord)
    # best_wing.calculate_mass()

    print(Fore.BLUE + f'Wing mass = {2*round(best_wing.m, 2)} [kg]')
    wn_wing, x_wing, U_wing = wing_vibrations()
    n_rivets_0 = best_wing.design_joint(b=b/2)
    n_rivets_1 = best_wing.design_joint(b=0)

    print(Fore.WHITE + '\n### Tail sizing started ###\n')
    hStabilizer, vStabilizer, tailPoleMass = size_tail()
