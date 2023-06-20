from dse.detailed.Structures.StructureClasses import Beam, Force, xflr_forces, TailVibes
from dse.detailed.Structures.rotor_sizing import y_transformation
from dse.detailed.Structures.material_properties import materials
from dse.detailed.Structures.constants import *
from dse.plotting import format_plot, save_plot
import vibration_toolbox as vtb
import matplotlib.pyplot as plt
from colorama import Fore
import pandas as pd
import numpy as np


def size_structure():
    # Rotors
    print(Fore.WHITE + "\n### Rotor blade sizing started ###\n")
    frontBlade, rearBlade, mr = size_rotor_blades()
    # mr = 344

    # Wings
    print(Fore.WHITE + "\n### Wing sizing started ###\n")
    wing, f_loading, moments = size_wing(wingDims['span'], wingDims['rootChord'], wingDims['taper'], mr, -1)
    wing.m_loading = moments
    wing_vibrations(wing)
    wing.design_joint(b=np.min(wing.y) / 2)
    wing.design_joint(b=0)

    # Tails
    print(Fore.WHITE + "\n### Tail sizing started ###\n")
    hStabilizer, vStabilizer, tailPoleMass = size_tail()

    return frontBlade, rearBlade, wing, hStabilizer, vStabilizer


def size_rotor_blades(overwrite=False):
    discretization = 5
    frontTwist = np.zeros(discretization * (np.size(rotorDims["frontBladeTwist"]) - 1))
    rearTwist = np.zeros(discretization * (np.size(rotorDims["rearBladeTwist"]) - 1))
    for i in range(np.size(rotorDims["frontBladeTwist"]) - 1):
        frontTwist[i * discretization : i * discretization + discretization] = np.linspace(
            rotorDims["frontBladeTwist"][i], rotorDims["frontBladeTwist"][i + 1], discretization
        )
        rearTwist[i * discretization : i * discretization + discretization] = np.linspace(
            rotorDims["rearBladeTwist"][i], rotorDims["rearBladeTwist"][i + 1], discretization
        )

    frontTwist = frontTwist[
        : np.size(frontTwist) - round(rotorDims["cutout"] * np.size(frontTwist))
    ]
    l_front = np.linspace(
        -rotorDims["radius"], -rotorDims["radius"] * rotorDims["cutout"], np.size(frontTwist)
    )

    rearTwist = rearTwist[: np.size(rearTwist) - round(rotorDims["cutout"] * np.size(rearTwist))]
    l_rear = np.linspace(
        -rotorDims["radius"], -rotorDims["radius"] * rotorDims["cutout"], np.size(frontTwist)
    )

    frontChord = np.linspace(
        rotorDims["chord"][0],
        rotorDims["chord"][
            np.size(rotorDims["chord"]) - round(rotorDims["cutout"] * np.size(rotorDims["chord"]))
        ],
        np.size(frontTwist),
    )

    rearChord = np.linspace(
        rotorDims["chord"][0],
        rotorDims["chord"][
            np.size(rotorDims["chord"]) - round(rotorDims["cutout"] * np.size(rotorDims["chord"]))
        ],
        np.size(rearTwist),
    )

    frontSect = rearSect = np.vstack((const["Airfoil"]["x"], const["Airfoil"]["z"]))

    frontSect = (
        np.ones((np.size(l_front), 2, 1))
        * frontSect
        * np.reshape(frontChord, (np.size(frontChord), 1, 1))
    )
    frontSect = y_transformation(frontTwist, frontSect)

    rearSect = (
        np.ones((np.size(l_rear), 2, 1))
        * rearSect
        * np.reshape(rearChord, (np.size(rearChord), 1, 1))
    )
    rearSect = y_transformation(rearTwist, rearSect)

    global Xac_rotor, Zac_rotor
    Xac_rotor = np.max(const["Airfoil"]["x"]) * frontChord[-1] / 4
    Zac_rotor = 0.077 * frontChord[0] / 20

    frontBlade = Beam(
        width=frontSect[:, 0].T,
        height=frontSect[:, 1].T,
        length=l_front,
        cross_section=frontSect,
        material="CFRCy",
        fixing_points=np.array([[Xac_rotor], [Zac_rotor]]) * np.ones(np.size(l_front)),
        name='front_blade'
    )

    rearBlade = Beam(
        width=rearSect[:, 0].T,
        height=rearSect[:, 1].T,
        length=l_rear,
        cross_section=rearSect,
        material="CFRCy",
        fixing_points=np.array([[Xac_rotor], [Zac_rotor]]) * np.ones(np.size(l_front)),
        name='rear_blade'
    )

    # Define the applied forces
    liftOffperBlade = (
        np.ones(np.shape(frontTwist))
        * 3071.90776835 / const["bladePerRotor"]
        / np.size(frontTwist)
    )
    liftOffDrag = (
        np.ones(np.shape(frontTwist)) * 250.676835 / (rotorDims['radius'] * (1 - rotorDims['cutout'])) / np.size(frontTwist)
    )
    application = np.ones(np.shape(frontTwist)) * np.array(
        [[Xac_rotor], [-rotorDims["radius"]], [Zac_rotor]]
    )
    application[1] = l_front

    liftOffForce_front = Force(
        magnitude=liftOffperBlade * np.array([[0], [0], [1]]),
        point_of_application=application,
    )

    liftOffDrag_front = Force(
        magnitude=liftOffDrag * np.array([[-1], [0], [0]]),
        point_of_application=application,
    )

    application2 = np.ones(np.shape(rearTwist)) * np.array(
        [[Xac_rotor], [-rotorDims["radius"]], [Zac_rotor]]
    )
    application2[1] = l_rear
    liftOffForce_rear = Force(
        magnitude=liftOffperBlade * np.array([[0], [0], [1]]),
        point_of_application=application2,
    )
    liftOffDrag_rear = Force(
        magnitude=liftOffDrag * np.array([[-1], [0], [0]]),
        point_of_application=application2,
    )

    diff = 100
    rotorMass = 250 / 9 / np.size(l_front)
    frontRotorMass = rotorMass
    rearRotorMass = rotorMass
    while diff > 0.01:
        frontBlade.unload()
        rearBlade.unload()

        rotationForce_front = frontRotorMass * const["rpm"] ** 2 * -l_front
        rotatingForce_front = Force(
            magnitude=np.vstack(
                (
                    np.zeros(np.size(rotationForce_front)),
                    rotationForce_front,
                    np.zeros(np.size(rotationForce_front)),
                )
            ),
            point_of_application=application,
        )

        rotationForce_rear = rearRotorMass * const["rpm"] ** 2 * -l_rear
        rotatingForce_rear = Force(
            magnitude=np.vstack(
                (
                    np.zeros(np.size(rotationForce_rear)),
                    rotationForce_rear,
                    np.zeros(np.size(rotationForce_rear)),
                )
            ),
            point_of_application=application,
        )

        frontBlade.add_loading(liftOffForce_front)
        frontBlade.add_loading(liftOffDrag_front)
        frontBlade.add_loading(rotatingForce_front)
        frontBlade.plot_internal_loading('front_blade')

        rearBlade.add_loading(liftOffForce_rear)
        rearBlade.add_loading(liftOffDrag_rear)
        rearBlade.add_loading(rotatingForce_rear)
        rearBlade.plot_internal_loading('rear_blade')

        frontBlade.InternalStress(0, 0, 0, x_scale=0.5, y_scale=0.5)
        frontBlade.m = np.sum(frontBlade.masses())

        rearBlade.InternalStress(0, 0, 0, x_scale=0.5, y_scale=0.5)
        rearBlade.m = np.sum(frontBlade.masses())

        diff = np.maximum(
            np.abs(rotorMass - frontBlade.m - 2) / (rotorMass - 2),
            np.abs(rotorMass - rearBlade.m - 2) / (rotorMass - 2),
        )
        frontRotorMass = np.hstack((np.sum(frontBlade.masses(), 0), np.array([0]))) + 2 * np.ones(
            np.size(l_front)
        ) / np.size(l_front)
        rearRotorMass = np.hstack((np.sum(rearBlade.masses(), 0), np.array([0]))) + 2 * np.ones(
            np.size(l_rear)
        ) / np.size(l_rear)

    def I(w, t):
        return w * t**3 + 2 * w * t * (0.073 - t) / 2 + t * (0.073 - 2 * t) ** 3 / 6

    required_I = 2.300534277610573e-06
    for i in np.linspace(0.001, frontChord[-1], 1000):
        if I(i, 0.001) > required_I:
            break

    print(
        f"A hollow square beam of thickness 1[mm] and width {round(i*1e3, 2)} [mm] is needed before the cutout"
    )
    rod_weight_1 = (
        (2 * i * 0.001 + 2 * 0.001 * (0.0840609 - 2 * 0.001))
        * 0.5
        * rotorDims["radius"]
        * materials["CFRCy"].rho
    )
    print(f"This beam weights {rod_weight_1} [kg]")

    print(
        f"Each front blade weights {np.round(frontBlade.m, 2) + 2} kg, including 2 kg of reinforcements"
    )
    print(
        f"Each rear blade weights {np.round(rearBlade.m, 2) + 2} kg, including 2 kg of reinforcements"
    )
    print(Fore.BLUE + f"Total rotor mass = {np.round(12 * (frontBlade.m + 2) + 12 * (rearBlade.m + 2), 2)} kg" + Fore.WHITE)

    if not overwrite:
        wn_rotor, x_rotor, U_rotor = rotor_vibrations(frontBlade)
        wn_rotor, x_rotor, U_rotor = rotor_vibrations(rearBlade)
    m_prop = 12 * (frontBlade.m + 2) + 12 * (rearBlade.m + 2) + rod_weight_1
    print(Fore.BLUE + f"The total mass of the propulsion subsystem is {round(m_prop, 2)} [kg]")
    return frontBlade, rearBlade, m_prop


def plot_mode_response(x, U):
    fig, ax = plt.subplots()
    ax.plot(x, U * 1e3, label=[f"Mode {i + 1}" for i in range(np.shape(U)[1])])
    ax.set_xlabel("Span [m]")
    ax.set_ylabel("Displacement [mm]")
    ax.legend()

    format_plot()
    save_plot('.', 'vibration_modes_wing')
    plt.show()


def equivalent_load(deflection, position, E, I):
    radius_of_curvature = position**2 / deflection
    curvature = 1 / radius_of_curvature

    # Assuming linearly elastic:
    M = curvature * E * I
    P = M / position
    return P


def rotor_vibrations(rotorBlade, reinforce=True, overwrite_I=None):
    L1 = np.abs(np.min(rotorBlade.y))
    E1 = materials["CFRCy"].E
    if np.shape(rotorBlade.Bi) == np.shape(rotorBlade.z[:-1]):
        I0 = np.mean(np.sum((rotorBlade.Bi * rotorBlade.z[:-1] ** 2), 0))
    else:
        I0 = np.mean(np.sum((rotorBlade.Bi * rotorBlade.z**2), 0))
    if overwrite_I is not None:
        I0 = overwrite_I
    A0 = np.mean(np.sum(rotorBlade.Bi, 0))

    print("\nVibration data:")
    print(f"Average moment of inertia of the airfoil = {I0}")

    bc = 2  # Clamped - free
    modes = 3

    parameters = np.array([E1, I0, materials["CFRCy"].rho, A0, L1])
    w, x, U = vtb.euler_beam_modes(n=modes, bctype=bc, beamparams=parameters)

    print(f"Original natural freq: {np.min(w)}")

    if reinforce:
        reinforcement_area = 0
        z_NA = np.sum(rotorBlade.Bi * rotorBlade.z[:-1], 0) / np.sum(rotorBlade.Bi, 0)
        z_max = np.max(rotorBlade.z, 0)
        arm = np.mean(z_max - z_NA)
        dA = 0.0001
        while np.min(w) < 25:
            A0 += dA
            I0 += dA * arm**2
            reinforcement_area += dA
            if reinforcement_area >= 1:
                raise ValueError("The added area is too large")

            parameters = np.array([E1, I0, materials["CFRCy"].rho, A0, L1])
            w, x, U = vtb.euler_beam_modes(n=modes, bctype=bc, beamparams=parameters)

        print(Fore.BLUE + f"Natural frequencies = {w} [rad/s]" + Fore.WHITE)
        print(f"Maximum deflection = {np.max(np.abs(U))} [m]")
        print(f"Required reinforcement area = {reinforcement_area}")
        print(
            f'Required reinforcement mass = {reinforcement_area*rotorDims["radius"]*(1-rotorDims["radius"])*materials["CFRCy"].rho} kg'
        )

        # Calculate equivalent load and additional stress
        max_deflection = U[np.where(np.abs(U) == np.max(np.abs(U)))]
        max_deflection_pos = x[np.where(np.abs(U) == np.max(np.abs(U)))[0]]
        P = equivalent_load(deflection=max_deflection, position=max_deflection_pos, E=E1, I=I0)

        print(f"Equivalent load due to vibrations in rotor blades is {P} [N]")

    plot_mode_response(x, U)
    return w, x, U


def wing_vibrations(wing, pars=None):
    # Define parameters
    if pars is None:
        L = -np.min(wing.y)
        I = np.mean(np.sum((wing.Bi * wing.z[:-1] ** 2), 0))
        A = np.mean(np.sum(wing.Bi, 0))
        parameters = np.array([materials["CFRPeek"].E, I, materials["CFRPeek"].rho, A, L])
    else:
        parameters = pars

    # Define analysis conditions
    bc = 3  # Clamped-pinned
    modes = 3

    # Perform analysis
    w, x, U = vtb.euler_beam_modes(n=modes, bctype=bc, beamparams=parameters)

    # Report results
    print(Fore.WHITE + "\nVibration data:")
    print(f"Average moment of inertia of the airfoil = {parameters[1]}")
    print(Fore.BLUE + f"Natural frequency of the wing = {w}")
    print(Fore.WHITE + f"Max deflection = {np.max(np.abs(U))} [m]")
    plot_mode_response(x, U)
    return w, x, U


def size_wing(span, chord_root, taper, rotor_mass=500, wing_model=None):
    mr = rotor_mass
    l = np.linspace(-span, 0, 100)
    chord_array = np.linspace(chord_root * taper, chord_root, np.size(l))

    # Define the geometry
    global Xac_wing, Zac_wing
    Xac_wing = np.max(const["Airfoil"]["x"]) * chord_array[-1] / 4
    Zac_wing = 0.077 * chord_array[-1]

    section = np.vstack((const["Airfoil"]["x"], const["Airfoil"]["z"])) * np.reshape(
        chord_array, (np.size(chord_array), 1, 1)
    )

    wing = Beam(
        width=np.array([const["Airfoil"]["x"]]).T * chord_array,
        height=np.array([const["Airfoil"]["z"]]).T * chord_array,
        length=l,
        cross_section=section,
        material="CFRPeek",
        fixing_points=np.array([[Xac_wing], [Zac_wing]]) * np.ones(np.size(l)),
        name='wing'
    )

    theta = np.arctan(const["fuselageHeight"] / span)

    # Define the forces during TO
    bracing_TO_mag = (
        const["g"] / np.sin(theta) * (1.1 * const["MTOM"] / 2 - (mr / 2 + const["engineMass"]))
    )
    bracing_TO = Force(
        magnitude=bracing_TO_mag * np.array([[0], [-np.cos(theta)], [-np.sin(theta)]]),
        point_of_application=np.array([[Xac_wing], [-span], [Zac_wing]]),
    )

    R_brace = bracing_TO_mag / (2 * np.pi * 0.001 * materials['CFRPeek'].compressive / 4.5)
    m_brace = materials['CFRPeek'].rho * span/np.cos(theta) * 2 * np.pi * R_brace * 0.001
    print(Fore.BLUE + f'The brace needs to have a radius of {R_brace} [m] and will weight {m_brace} [kg]')
    engine_and_rotor_weight = Force(
        magnitude=np.array([[0], [0], [-(mr / 2 + const["engineMass"]) * const["g"]]]),
        point_of_application=np.array([[Xac_wing], [-span], [Zac_wing]]),
    )
    liftOffLoad = Force(
        magnitude=np.array([[0], [0], [1.1 * const["MTOM"] * const["g"] / 2]]),
        point_of_application=np.array([[Xac_wing], [-span], [Zac_wing]]),
    )

    # Define loads during cruise
    if wing_model is None:
        aerodynamic_forces = xflr_forces("Test_xflr5_file.csv", const["q"], float(span))
    else:
        if wing_model == -1:
            cl, cd, y_L = xflr_forces("freek.csv", const["q"], float(span), adrian=wing_model)
        else:
            cl, cd, y_L = xflr_forces("wings.csv", const["q"], float(span), adrian=wing_model)

        dy = np.abs(wing.y[:-1] - wing.y[1:])
        c_arr = (chord_array[:-1] + chord_array[1:]) / 2

        Cl = np.zeros(np.size(dy))
        Cd = np.zeros(np.size(dy))
        app = np.zeros(np.size(dy))
        for i in range(np.size(dy)):
            indx = (np.abs(wing.y[i] + y_L)).argmin()
            Cl[i] = cl[indx]
            Cd[i] = cd[indx]
            app[i] = wing.y[i]

        lift = const["load_factor"] * Cl * const["q"] * c_arr * dy
        drag = const["load_factor"] * Cd * const["q"] * c_arr * dy

        aerodynamic_forces = Force(
            magnitude=np.vstack((-drag, np.zeros(np.size(lift)), lift)),
            point_of_application=np.vstack(
                (Xac_wing * np.ones(np.size(lift)), app, Zac_wing * np.ones(np.size(lift)))
            ),
        )

    cruiseThrust = Force(
        magnitude=np.array([[-np.sum(aerodynamic_forces.F[0])], [0], [0]]),
        point_of_application=np.array([[Xac_wing], [-span], [Zac_wing]]),
    )

    liftMoment = -np.dot(aerodynamic_forces.F[2], aerodynamic_forces.application[1])
    bracing_cr_mag = (
        1
        / (span * np.sin(theta))
        * (liftMoment - span * const["g"] * (mr / 2 + const["engineMass"]))
    )
    bracing_cruise = Force(
        magnitude=bracing_cr_mag * np.array([[0], [-np.cos(theta)], [-np.sin(theta)]]),
        point_of_application=np.array([[Xac_wing], [-span], [Zac_wing]]),
    )


    # Apply loads and size
    wing.unload()
    wing.add_loading(liftOffLoad)
    wing.add_loading(engine_and_rotor_weight)
    wing.add_loading(bracing_TO)
    wing.plot_internal_loading('wing_TO')
    wing.InternalStress(0, 0, 0)
    f_loading_TO = wing.f_loading
    thickness = wing.t
    B1 = wing.Bi
    stress_TO = wing.sigma

    # Apply loads and size
    wing.add_loading(engine_and_rotor_weight)
    wing.add_loading(cruiseThrust)
    wing.add_loading(bracing_cruise)
    wing.add_loading(aerodynamic_forces)
    wing.plot_internal_loading(structure='wing_cruise')
    wing.InternalStress(0, 0, 0, y_scale=0.75)
    f_loading_Cr = wing.f_loading
    moments = wing.m_loading
    thickness2 = wing.t
    B2 = wing.Bi
    stress_cr = wing.sigma

    # Choose the most critical structure
    wing.t = np.maximum(thickness, thickness2)
    wing.Bi = np.maximum(B1, B2)
    wing.sigma = np.minimum(stress_TO, stress_cr)
    wing.sigma[np.where(stress_cr > 0)] = np.maximum(stress_cr[np.where(stress_cr > 0)], np.abs(stress_TO[np.where(stress_cr > 0)]))
    f_loading_abs = np.maximum(np.abs(f_loading_TO), np.abs(f_loading_Cr))

    print(f'Buckling in cruise:')
    wing.calculate_mass()
    print(Fore.BLUE + f"Mass of each wing = {np.round(wing.m, 2)}kg" + Fore.WHITE)
    return wing, f_loading_abs, moments


def size_tail():
    # span, 14.23/2  chord 2.7106, tip = 1.0834 NACA0012
    # Assumptions
    m = 25

    ### Horizontal stabilizer ###
    # Define the geometry
    tailChord = np.linspace(hTailDims["tipChord"], hTailDims["rootChord"], m)
    NACA0012 = pd.read_csv(
        ".\\dse\\detailed\\Structures\\NACA 0012.dat", delimiter="\s+", dtype=float, skiprows=1, names=["x", "z"]
    )

    l = np.linspace(-hTailDims["span"], 0, m)
    x = np.reshape(NACA0012["x"], (len(NACA0012["x"]), 1))
    z = np.reshape(NACA0012["z"], (len(NACA0012["z"]), 1))
    section = np.vstack((NACA0012["x"], NACA0012["z"])) * np.reshape(tailChord, (m, 1, 1))

    Xac = np.max(NACA0012["x"]) * tailChord[0] / 4
    Zac = 0.077 * tailChord[0]

    hStabilizer = Beam(
        width=x * tailChord,
        height=z * tailChord,
        length=l,
        cross_section=section,
        material="CFRPeek",
        fixing_points=np.array([[Xac], [Zac]]) * np.ones(m),
        name='hStabilizer'
    )

    # Define the loads
    lift = Force(  # Assuming uniform distribution for now, will change later
        magnitude=(const["tailWingLift"] * const["MTOM"] * const["g"] + const["extraForceTail"])
        / m
        * np.vstack((np.zeros((2, m)), np.ones(m))),
        point_of_application=np.vstack((Xac * np.ones(m), l, Zac * np.ones(m))),
    )

    drag = Force(  # Assuming uniform distribution for now, will change later
        magnitude=(const["tailWingLift"] * const["MTOM"] * const["g"] + const["extraForceTail"])
        / (m * const["liftToDragTail"])
        * np.vstack((-1 * np.ones(m), np.zeros((2, m)))),
        point_of_application=np.vstack((Xac * np.ones(m), l, Zac * np.ones(m))),
    )

    # Apply loads and size
    hStabilizer.add_loading(lift)
    hStabilizer.add_loading(drag)
    hStabilizer.plot_internal_loading('hStabilizer')
    hStabilizer.InternalStress(0, 0, 0)
    hStabilizer.calculate_mass()
    print(f"Horizontal stabilizer's mass is {2*np.round(hStabilizer.m)} [kg]")

    ### Vertical stabilizer ###
    # Define geometry
    vTailChord = np.linspace(vTailDims["tipChord"], vTailDims["rootChord"], m)

    lv = np.linspace(-vTailDims["span"], 0, m)
    xv = np.reshape(NACA0012["x"], (len(NACA0012["x"]), 1))
    zv = np.reshape(NACA0012["z"], (len(NACA0012["x"]), 1))
    section = np.vstack((NACA0012["x"], NACA0012["z"])) * np.reshape(vTailChord, (m, 1, 1))

    Xac = np.max(NACA0012["x"]) * vTailChord[0] / 4
    Zac = 0.077 * vTailChord[0]

    vStabilizer = Beam(
        width=xv * vTailChord,
        height=zv * vTailChord,
        length=lv,
        cross_section=section,
        material="AL_light",
        fixing_points=np.array([[Xac], [Zac]]) * np.ones(m),
        name='vStabilizer'
    )

    # Define loads
    restoring = Force(
        magnitude=const["extraForceTail"]
        / m
        * np.vstack((-np.ones(m) / const["liftToDragTail"], np.zeros(m), np.ones(m))),
        point_of_application=np.vstack(
            (
                Xac * np.ones(m),
                lv,
                Zac * np.ones(m),
            )
        ),
    )

    # Apply loads and size
    vStabilizer.add_loading(restoring)
    vStabilizer.plot_internal_loading('vStabilizer')
    vStabilizer.InternalStress(0, 0, 0)
    vStabilizer.calculate_mass()
    print(f"Vertical stabilizer's mass is {np.round(vStabilizer.m)} [kg]")

    ### Tail pole ###
    boomMoments = vStabilizer.m_loading[-1] + const["tailPoleLength"] * (
        np.array([[-hStabilizer.f_loading[-1][2][0]], [0], [0]])
        + np.array([[0], [0], [-vStabilizer.f_loading[-1][2][0]]])
        + np.array([[const["g"] * (hStabilizer.m + vStabilizer.m)], [0], [0]])
    )

    s_max = materials["CFRP"].compressive
    t_min = 0.001
    Mx = boomMoments[0]
    Mz = boomMoments[2]
    Fy = hStabilizer.f_loading[-1][0][0] + vStabilizer.f_loading[-1][1][0]
    theta = np.arctan(Mx / Mz)
    D = (
        Fy
        + np.sqrt(
            Fy**2 + 8 * s_max * np.pi * t_min * np.abs(Mz * np.cos(theta) + Mx * np.sin(theta))
        )
    ) / (2 * s_max * np.pi * t_min)

    tailPoleMass = np.pi * D * 0.001 * const["tailPoleLength"] * materials["CFRPeek"].rho

    print(f"Tail pole mass = {tailPoleMass} [kg]")
    print(
        Fore.BLUE
        + f"The total tail group mass is {2*hStabilizer.m + vStabilizer.m + tailPoleMass} [kg]"
    )
    # l = 15 - 2.1439
    # t = 0.001
    # r = tailPoleMass / (2 * np.pi * t * l * materials['CFRCy'].rho82)
    # tail_mass = hStabilizer.m + vStabilizer.m + tailPoleMass
    # Th = TailVibes(E=materials['CFRCy'].E, density=materials['CFRCy'].rho, radius=r, length=l, thickness=t,
    #               tail_mass=tail_mass, surface_area=16.955)
    #
    # Th.simsetting()
    # Th.sysparam()
    # Th.userinput(cl=1)
    # Th.syssim()
    # Th.results()
    #
    # Tv = TailVibes(E=materials['CFRCy'].E, density=materials['CFRCy'].rho, radius=r, length=l, thickness=t,
    #               tail_mass=tail_mass, surface_area=11.333)
    #
    # Tv.simsetting()
    # Tv.sysparam()
    # Tv.userinput(cl=1)
    # Tv.syssim()
    # # T.results()
    # nth = int((1 / Th.dt) / 1e4)
    # plt.figure(figsize=(9, 2.5))
    # plt.plot(Th.t[::nth], Th.x[::nth], label="Displacement due to elevator deflection")
    # plt.plot(Tv.t[::nth], Tv.x[::nth], label="Displacement due to rudder deflection", linestyle="dashed")
    # plt.legend()
    # # plt.plot(self.t[::nth], self.v[::nth], label="Velocity")
    # # plt.plot(self.t[::nth], self.F_u[::nth]/self.keq, label="Force")
    # # plt.axhline(self.avg)
    # # plt.legend()
    # plt.xlabel('Time [s]')
    # plt.ylabel('Displacement [m]')
    # format_plot()
    # save_plot('.', 'Tail vibrations')
    # plt.show()
    return hStabilizer, vStabilizer, tailPoleMass


if __name__ == "__main__":
    structures = size_structure()
