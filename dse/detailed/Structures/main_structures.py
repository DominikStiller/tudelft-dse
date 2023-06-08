from StructureClasses import Beam, Force, xflr_forces
from rotor_sizing import y_transformation
from material_properties import materials
import vibration_toolbox as vtb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


def size_rotor_blades():
    print(f'\n### Rotor blade sizing started ###\n')
    # Define the blade
    bladeTwist = np.flip(np.radians(np.array([20.01066435, 15.27067334, 12.97291975, 11.55038814, 10.55750286,
                                              9.81326113, 9.22850394,  8.75360237,  8.35851135,  8.02392466,
                                              7.73691493,  7.48862104, 7.27297283,  7.08600653,  6.92558553,
                                              6.79150814,  6.6861861,   6.61653976, 6.59931193, 6.67922986])))
    rootBladeChord = 0.693
    tipBladeChord = 0.347
    cutout = 0.5
    discretization = 20

    twist = np.zeros(discretization * (np.size(bladeTwist)-1))
    for i in range(np.size(bladeTwist)-1):
        twist[i*discretization:i*discretization+discretization] = np.linspace(bladeTwist[i], bladeTwist[i+1], discretization)


    twist = twist[:np.size(twist)-round(cutout*np.size(twist))]
    l = np.linspace(-R, -R * cutout, np.size(twist))

    bladeChord = np.linspace(tipBladeChord, rootBladeChord, np.size(twist))
    sect = np.vstack((Airfoil['x'], Airfoil['z']))
    sect = np.ones((np.size(twist), 2, np.size(Airfoil['x']))) * sect
    sect = y_transformation(twist, sect) * np.reshape(bladeChord, (np.size(bladeChord), 1, 1))

    Xac = np.max(Airfoil['x']) * rootBladeChord / 4
    Zac = 0.077 * tipBladeChord / 20

    blade = Beam(
        width=sect[:, 0].T,
        height=sect[:, 1].T,
        length=l,
        cross_section=sect,
        material='CFRP',
        fixing_points=np.array([[Xac], [Zac]]) * np.ones(np.size(l))
    )

    # Define the applied forces
    liftOffperBlade = np.ones(np.shape(twist)) * 1.1 * MTOM * g / 24 / np.size(twist)
    application = np.ones(np.shape(twist)) * np.array([[Xac], [-R], [Zac]])
    application[1] = l

    liftOffForce = Force(
        magnitude=liftOffperBlade * np.array(
            [
                [-1 / 54.178],
                [0],
                [1]
            ]
        ),
        point_of_application=application
    )

    diff = 100
    rotorMass = 250 / 9 / np.size(l)
    while diff > 0.01:
        blade.unload()

        rotationForce = rotorMass * rpm**2 * l
        rotatingForce = Force(
            magnitude=np.vstack((np.zeros(np.shape(rotationForce)), rotationForce, np.zeros(np.shape(rotationForce)))),
            point_of_application=application
        )

        blade.add_loading(liftOffForce)
        blade.add_loading(rotatingForce)
        blade.plot_internal_loading()

        blade.InternalStress(0, 0, 0)
        blade.calculate_mass()

        diff = np.abs(rotorMass - blade.m - 10) / (rotorMass - 10)
        rotorMass = np.hstack((np.sum(blade.masses(), 0), np.array([0]))) + 10 * np.ones(np.size(l)) / np.size(l)

    boomMoments = blade.m_loading[-1] + R*cutout * np.array(
        [
            [blade.f_loading[-1][2][0]],
            [0],
            [blade.f_loading[-1][0][0]]
        ]
    )

    Mx = boomMoments[0]
    Mz = boomMoments[2]
    Fy = blade.f_loading[-1][1][0]
    maxStressPos = np.arctan(Mx / Mz)
    rodMat = materials['Titanium Alloys']
    D = (Fy + np.sqrt(
        Fy ** 2 + 8 * rodMat.compressive / 1.5 * np.pi * 0.001 * (Mz * np.cos(maxStressPos) + Mx * np.sin(maxStressPos)))) / (
                2 * rodMat.compressive / 1.5 * np.pi * 0.001)

    rod_weight = np.pi * D * 0.001 * R * rodMat.rho
    print(f'Minimum rod diameter = {1000*D} [mm], with a weight of {rod_weight} [kg]')

    m_r = rotorMass * 24

    print(f'Each blade weights {np.round(blade.m, 2) + 10} kg, including kg of reinforcements')
    print(f'Total rotor mass = {np.round(24 * (blade.m + 10), 2)} kg')

    return blade, np.sum(m_r), D


def plot_rotor_vibrations():
    # Parameters of the rod
    cutout = 0.5
    L1 = R
    E1 = materials['CFRP'].E
    I0 = np.mean(np.sum((rotorBlade.Bi * rotorBlade.z[:-1]**2), 0))
    A0 = np.mean(np.sum(rotorBlade.Bi, 0))

    print('\nVibration data:')
    print(f'Average moment of inertia of the airfoil = {I0}')

    bc = 3  # Clamped-pinned
    modes = 3

    parameters = np.array([E1, I0, materials['CFRP'].rho, A0, L1])
    w, x, U = vtb.euler_beam_modes(n=modes, bctype=bc, beamparams=parameters)

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

        parameters = np.array([E1, I0, materials['CFRP'].rho, A0, L1])
        w, x, U = vtb.euler_beam_modes(n=modes, bctype=bc, beamparams=parameters)

    print(f'Natural frequencies = {w} [Hz]')
    print(f'Maximum deflection = {np.max(U)} [m]')
    print(f'Required reinforcement area = {reinforcement_area}')
    print(f'Required reinforcement mass = {reinforcement_area*R*(1-cutout)*materials["CFRP"].rho} kg')

    # Calculate equivalent load and additional stress

    fig, ax = plt.subplots()
    ax.plot(x, U*1e3, label=[f'Mode {i+1}' for i in range(np.shape(U)[1])])
    ax.set_xlabel('Span [m]')
    ax.set_ylabel('Displacement [mm]')

    ax.grid(True)
    ax.legend()

    fig.tight_layout()
    plt.show()


def size_wing(chord_array, span):
    # Define the geometry
    Xac = np.max(Airfoil['x']) * chord_array[-1] / 4
    Zac = 0.077 * chord_array[-1]

    l = np.linspace(-span, 0, 100)

    section = np.vstack((Airfoil["x"], Airfoil["z"])) * np.reshape(np.vstack(chord_array, chord_array), (np.size(chord_array), 2, 1))

    wing = Beam(
        width=Airfoil["x"] * chord_array,
        height=Airfoil["z"] * chord_array,
        length=l,
        cross_section=section,
        material='CFRP',
        fixing_points=np.array([[Xac], [Zac]]) * np.ones(np.size(l))
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
                [Xac],
                [-span],
                [Zac]
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
                [Xac],
                [-span],
                [Zac]
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
                [Xac],
                [-span],
                [Zac]
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
                [Xac],
                [-span],
                [Zac]
            ]
        )
    )

    # Apply loads and size
    wing.add_loading(liftOffLoad)
    wing.add_loading(engine_and_rotor_weight)
    wing.add_loading(bracing_TO)
    wing.plot_internal_loading()
    wing.InternalStress(0, 0, 0)
    thickness = wing.t
    B1 = wing.Bi

    # Define loads during cruise
    wing.unload()
    aerodynamic_forces = xflr_forces('Test_xflr5_file.csv', q, span)

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
                [Xac],
                [-span],
                [Zac]
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
    thickness2 = wing.t
    B2 = wing.Bi

    # Choose the most critical structure
    wing.t = np.maximum(thickness, thickness2)
    wing.Bi = np.maximum(B1, B2)

    wing.calculate_mass()
    print(f'Mass of each wing = {np.round(wing.m, 2)}kg')
    return wing


def size_tail():
    # Assumptions
    m = 25
    tail_to_wing_lift = 0.1
    cl_t = 1
    AR_tail = 10
    lift_to_drag_tail = 15
    extra_force = 500  # For control
    tail_taper = np.linspace(1, 0.7, m)
    vTailTaper = tail_taper

    ### Horizontal stabilizer ###
    # Define the geometry
    tailChord = tail_to_wing_lift * MTOM * g / (AR_tail * cl_t * q)
    tailSpan = tailChord * AR_tail / 2

    l = np.linspace(-tailSpan, 0, m)
    x = np.reshape(Airfoil['x'] * tailChord, (np.size(Airfoil['x']), 1))
    z = np.reshape(Airfoil['z'] * tailChord, (np.size(Airfoil['z']), 1))
    section = np.vstack((x.T, z.T)) * np.ones((m, 2, 1)) * np.reshape(tail_taper, (m, 1, 1))

    Xac = np.max(Airfoil['x']) * tailChord / 4
    Zac = 0.077 * tailChord

    hStabilizer = Beam(
        width=x * tail_taper,
        height=z * tail_taper,
        length=l,
        cross_section=section,
        material='CFRP',
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
    print(f"Horizontal stabilizer's mass is {np.round(hStabilizer.m)} [kg]")

    ### Vertical stabilizer ###
    # Define geometry

    vTailSpan = extra_force / (cl_t * q * AR_tail)
    vTailChord = vTailSpan / vTailSpan


    l = np.linspace(-vTailSpan, 0, m)
    x = np.reshape(Airfoil['x'] * vTailChord, (np.size(Airfoil['x']), 1))
    z = np.reshape(Airfoil['z'] * vTailChord, (np.size(Airfoil['z']), 1))
    section = np.vstack((x.T, z.T)) * np.ones((m, 2, 1)) * np.reshape(tail_taper, (m, 1, 1))

    Xac = np.max(Airfoil['x']) * vTailChord / 4
    Zac = 0.077 * vTailChord

    vStabilizer = Beam(
        width=x * vTailTaper,
        height=z * vTailTaper,
        length=l,
        cross_section=section,
        material='CFRP',
        fixing_points=np.array([[Xac], [Zac]]) * np.ones(m)
    )

    # Define loads
    restoring = Force(
        magnitude=extra_force / m * np.vstack((-np.ones(m)/lift_to_drag_tail, np.zeros(m), np.ones(m))),
        point_of_application=np.vstack(
            ((Xac * np.min(vTailTaper) + (1-np.min(vTailTaper))*vTailChord) * np.ones(m),
             l,
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

    tailPoleMass = np.pi * D * 0.001 * R * materials['CFRP'].rho

    print(f'Tail pole mass = {tailPoleMass} [kg]')
    print(f'Tail group mass = {hStabilizer.m + vStabilizer.m + tailPoleMass + margin} [kg], '
          f'including {margin} [kg] of margin')
    return hStabilizer, vStabilizer, tailPoleMass


def size_body(fuselage_height=1.67, cabin_length=2, full_length=6.15):
    r = fuselage_height/2
    aft_cone_length = full_length - cabin_length - rootChord  # [m], assumed
    mat = materials['CFRP']
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

    # Start sizing
    # Size the rotor blades for their loads
    rotorBlade, mr, D = size_rotor_blades()
    rotorBlade.overall_inertia()

    # Evaluate vibrational response
    plot_rotor_vibrations()

    span = np.array([4500, 4000, 3500, 3000])
    rootChord = np.array([4, 4.3333, 4.9481, 5.78])

    # wing = size_wing()
    # wing.overall_inertia()
    # nf = wing.natural_frequency()
    # print(f'Lowest natural frequency of the wing = {np.min(np.abs(nf[:, 1][np.nonzero(nf[:, 1])]))}')
    #
    #
    # hStabilizer, vStabilizer, tailPoleMass = size_tail()
    # hStabilizer.overall_inertia()
    # vStabilizer.overall_inertia()
    #
    # bodyMass, Ix_f, Iy_f, Iz_f = size_body()
    #
    # # Calculate Ixx of the AC
    #
    # # Wing contribution - Assume x and z coordinates of wing's cg coincide with aircraft's
    # Ix = 2 * (wing.Ix + wing.m * wing.ycg**2)
    # Iy = 2 * (wing.Iy + wing.m * (1.67/2)**2)
    # Iz = 2 * (wing.Iz + wing.m * wing.ycg**2)
    #
    # # Engine contribution - Assume x and z coordinates of engines' cg coincide with aircraft's
    # Ix += 2 * (2*m_e * b**2)
    # Iz += 2 * (2*m_e * b**2)
    #
    # # Rotor blade contributions - Assume x and z coordinates of engines' cg coincide with aircraft's
    # Ix += 2 * (12 * rotorBlade.Ix + rotorBlade.m * b**2)
    # Iy += 24 * rotorBlade.Iy
    # Iz += 2 * (12 * rotorBlade.Iz + rotorBlade.m + b**2)
    #
    # # Fuselage contribution
    # Ix += Ix_f
    # Iy += Iy_f
    # Iz += Iz_f
    #
    # # Tail contribution
    # Ix += hStabilizer.Ix + vStabilizer.Ix
    # Iy += hStabilizer.Iy + vStabilizer.Iy + (hStabilizer.m + vStabilizer.m) * (tailPoleMass + rootChord)**2
    # Iz += hStabilizer.Iz + vStabilizer.Iz + (hStabilizer.m + vStabilizer.m) * (tailPoleMass + rootChord)**2
    #
    # print(f'Overall Ix, Iy, Iz = {Ix, Iy, Iz} [kgm^2]')
