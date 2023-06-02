from StructureClasses import Beam, Force, xflr_forces
from rotor_sizing import y_transformation
from material_properties import materials
import pandas as pd
import numpy as np



def size_rotor_blades():
    # Define the blade
    bladeTwist = np.flip(np.radians(np.array([5.72821286, 5.29796352, 4.97731679, 4.74279731, 4.58226976,
                                              4.49233726, 4.47898847, 4.56296873, 4.79802654, 5.34554722])))

    l = np.linspace(-R, -R / 3, np.size(bladeTwist))

    sect = np.vstack((Airfoil['x'], Airfoil['z']))
    sect = np.ones((np.size(bladeTwist), 2, np.size(Airfoil['x']))) * sect
    sect = y_transformation(bladeTwist, sect) * R / 20

    Xac = np.max(Airfoil['x']) * R / 20 / 4
    Zac = 0.077 * R / 20

    blade = Beam(
        width=sect[:, 0].T,
        height=sect[:, 1].T,
        length=l,
        cross_section=sect,
        material=materials['CFRP'],
        fixing_points=np.array([[Xac], [Zac]]) * np.ones(np.size(l))
    )

    # Define the applied forces
    liftOffperBlade = np.ones(np.shape(bladeTwist)) * 1.1 * MTOM * g / 24 / np.size(bladeTwist)
    application = np.ones(np.shape(bladeTwist)) * np.array([[1.603 / 3.35 * R / 20], [-R], [0.1742 / 3.35 * R / 20]])
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
    rotorMass = 250 / 9
    while diff > 0.01:
        blade.unload()

        rotationForce = rotorMass * rpm * l
        rotatingForce = Force(
            magnitude=np.vstack((np.zeros(np.shape(rotationForce)), rotationForce, np.zeros(np.shape(rotationForce)))),
            point_of_application=application
        )

        blade.add_loading(liftOffForce)
        blade.add_loading(rotatingForce)
        blade.plot_internal_loading()

        blade.InternalStress(0, 0, 0)
        blade.calculate_mass()
        print(blade.m)

        diff = np.abs(rotorMass - blade.m - 10) / (rotorMass - 10)
        rotorMass = blade.m + 10

    boomMoments = blade.m_loading[-1] + R/3 * np.array(
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
    D = (Fy + np.sqrt(
        Fy ** 2 + 8 * 4e8 / 1.5 * np.pi * 0.001 * (Mz * np.cos(maxStressPos) + Mx * np.sin(maxStressPos)))) / (
                2 * 4e8 / 1.5 * np.pi * 0.001)

    rod_weight = np.pi * D * 0.001 * R * materials['CFRP'].rho

    print(f'Each blade weights {np.round(blade.m, 2) + 10} kg, including {np.round(rod_weight, 2)} kg of rod and '
          f'{np.round(10 - rod_weight, 2)} kg for other reinforcements')
    print(f'Total rotor mass = {np.round(24 * (blade.m + 10), 2)} kg')

    m_r = (rotorMass + 10) * 24
    return blade, m_r


def size_wing():
    # Define the geometry
    Xac = np.max(Airfoil['x']) * rootChord / 4
    Zac = 0.077 * rootChord
    l = np.linspace(-b, 0, 100)

    wing = Beam(
        width=Airfoil["x"].to_numpy() * rootChord,
        height=Airfoil["z"].to_numpy() * rootChord,
        length=l,
        cross_section="constant",
        material=materials['CFRP'],
        fixing_points=np.array([[Xac], [Zac]])
    )
    theta = np.arctan(fuselage_height / b)

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
                [-b],
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
                [-b],
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
                [-b],
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
                [-b],
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
    aerodynamic_forces = xflr_forces('Test_xflr5_file.csv', q, b)

    liftMoment = -np.dot(aerodynamic_forces.F[2], aerodynamic_forces.application[1])
    bracing_cr_mag = 1 / (b * np.sin(theta)) * (liftMoment - b * g * (mr/2 + m_e))
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
                [-b],
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
    cabin_length = 1.73
    x_ac_w = cabin_length + rootChord / 4
    x_cg = x_ac_w + 1
    x_ac_t = x_cg + 10
    lengthTailPole = x_ac_t - cabin_length - rootChord

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
        material=materials['CFRP'],
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
        material=materials['CFRP'],
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


if __name__ == '__main__':
    # Constants
    R = 10.4
    MTOM = 3000
    g = 3.71
    rpm = 183 * 2 * np.pi / 60
    q = 0.5 * 0.01 * 112 ** 2
    rootChord = 3.35
    b = 16.8
    fuselage_height = 2
    m_e = 50

    Airfoil = pd.read_csv("S1223.dat", delimiter="\s+", dtype=float, skiprows=1, names=["x", "z"])

    rotorBlade, mr = size_rotor_blades()
    wing = size_wing()
    hStabilizer, vStabilizer, tailPoleMass = size_tail()
