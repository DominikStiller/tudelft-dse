from material_properties import  materials
from loading import Beam, Force
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def y_transformation(theta, arr):
    if len(np.shape(arr)) == 2:
        av = np.mean(arr, 1)
        T = np.matrix([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
        rotated = np.asarray(T * arr)
        rotated += np.reshape(av-np.mean(rotated, 1), (2, 1))
        return rotated
    elif len(np.shape(arr)) == 3:
        dummy = np.zeros(np.shape(arr))
        for i in range(np.shape(arr)[0]):
            av = np.mean(arr[i], 1)
            T = np.matrix([[np.cos(theta[i]), np.sin(theta[i])], [-np.sin(theta[i]), np.cos(theta[i])]])
            dummy[i] = T * arr[i]
            dummy[i] += np.reshape(av - np.mean(dummy[i], 1), (2, 1))
        return np.asarray(dummy)
    else:
        raise TypeError('The array needs to be either 2D or 3D')


if __name__ == '__main__':
    # Twist root-to-tip
    bladeTwist = np.flip(np.radians(np.array([5.72821286, 5.29796352,  4.97731679,  4.74279731,  4.58226976,
                                              4.49233726,  4.47898847, 4.56296873,  4.79802654, 5.34554722])))

    bladeTwistv = np.reshape(bladeTwist, (np.size(bladeTwist), 1, 1))
    twistCorrection = np.hstack((np.cos(bladeTwistv), np.sin(bladeTwistv)))

    R = 10.4
    Airfoil = pd.read_csv("S1223.dat", delimiter="\s+", dtype=float, skiprows=1, names=["x", "z"])
    l = np.linspace(-R, -R/3, np.size(bladeTwist))

    sect = np.vstack((Airfoil['x'], Airfoil['z']))
    sect = np.ones((np.size(bladeTwist), 2, np.size(Airfoil['x']))) * sect
    sect = y_transformation(bladeTwist, sect) * R/20

    plt.plot(sect[0][0], sect[0][1], label='root')
    plt.plot(sect[6][0], sect[6][1], label='middle')
    plt.legend()
    plt.gca().set_aspect('equal')
    plt.show()

    Xac = np.max(Airfoil['x']) * R/20 / 4
    Zac = 0.077 * R/20

    blade = Beam(
        width=sect[:, 0].T,
        height=sect[:, 1].T,
        length=l,
        cross_section=sect,
        material=materials['CFRP'],
        fixing_points=np.array([[Xac], [Zac]])*np.ones(np.size(l))
    )


    MTOM = 3000
    g = 3.71
    liftOffperBlade = np.ones(np.shape(bladeTwist)) * 1.1 * MTOM * g / 24 / np.size(bladeTwist)
    application = np.ones(np.shape(bladeTwist)) * np.array([[1.603/3.35*R/20], [-R], [0.1742/3.35*R/20]])
    application[1] = l

    liftOffForce = Force(
        magnitude=liftOffperBlade * np.array(
            [
                [-1/54.178],
                [0],
                [1]
            ]
        ),
        point_of_application=application
    )

    rpm = 183 * 2 * np.pi/ 60
    massGuess = 250/9

    print(massGuess)
    diff = 100
    while diff > 0.01:
        blade.unload()

        rotationForce = massGuess * rpm * l
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

        diff = np.abs(massGuess - blade.m - 10) / (massGuess-10)
        massGuess = blade.m + 10


    boomMoments = blade.m_loading[-1] + np.array(
        [
            [blade.f_loading[-1][0][0]],
            [0],
            [blade.f_loading[-1][2][0]]
        ]
    )

    Mx = boomMoments[0]
    Mz = boomMoments[2]
    Fy = blade.f_loading[-1][1][0]
    maxStressPos = np.arctan(Mx/Mz)
    D = (Fy + np.sqrt(
        Fy ** 2 + 8 * 4e8 / 1.5 * np.pi * 0.001 * (Mz * np.cos(maxStressPos) + Mx * np.sin(maxStressPos)))) / (
                    2 * 4e8 / 1.5 * np.pi * 0.001)

    rod_weight = np.pi*D*0.001*R*materials['CFRP'].rho

    print(f'Each blade weights {np.round(blade.m, 2)+10} kg, including {np.round(rod_weight, 2)} kg of rod and '
          f'{np.round(10-rod_weight, 2)} kg for other reinforcements')
    print(f'Total rotor mass = {np.round(24*(blade.m + 10), 2)} kg')
