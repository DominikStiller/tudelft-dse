import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit

from .constants import const
from .cruise_sizing import max_tipSpeed
from .power_sizing import power


def RadiusMassElementMomentum(
    M,
    N_rotors,
    N_blades,
    coaxial,
    V_tip,
    print_results=False,
    changeConstants=None,
    ExtremeValue=False,
):
    if changeConstants is not None:
        const[changeConstants[0]] = changeConstants[1]

    # Rename constants
    gm = const["gravityMars"]
    rho = const["airDensity"]
    V_m = const["soundSpeed"]

    if N_rotors <= 0:
        return "N_rotors has to be greater than zero.", 0, 0, 0
    elif N_blades <= 0:
        return "N_blades has to be greater than zero.", 0, 0, 0

    b = N_blades
    v_tip = V_tip
    T_min = M * gm / N_rotors  # Thrust to be produced per rotor

    # Initial guesses
    T = 0
    R = 1

    if ExtremeValue is True:
        c = R / 10  # Update the chord
        pow = power(T, R)

    while T < T_min:
        R += 0.1  # Increase the radius
        # Rename constants
        gm = const["gravityMars"]
        rho = const["airDensity"]
        V_m = const["soundSpeed"]

        # if N_rotors <= 0:
        #     return "N_rotors has to be greater than zero.",0,0,0
        # elif N_blades <= 0:
        #     return "N_blades has to be greater than zero.",0,0,0

        b = N_blades
        v_tip = V_tip

        c = R / 20  # Update the chord

        x0 = c / 2  # Distance from the blade's base to the centre of rotation
        omega = v_tip / R * 0.92  # Angular velocity

        n_elements = 15  # Split up blade into set of elements

        takeoffAOA = np.radians(6)

        alpha_data = [
            -5.5,
            -5.25,
            -5,
            -4.75,
            -4.5,
            -4.25,
            -4,
            -3.75,
            -3.5,
            -3.25,
            -3,
            -2.75,
            -2.5,
            -2.25,
            -2,
            -1.75,
            -1.5,
            -1.25,
            -1,
            -0.75,
            -0.5,
            -0.25,
            0,
            0.25,
            0.5,
            0.75,
            1,
            1.25,
            1.5,
            1.75,
            2,
            2.25,
            2.5,
            2.75,
            3,
            3.25,
            3.5,
            3.75,
            4,
            4.25,
            4.5,
            4.75,
            5,
            5.25,
            5.5,
            5.75,
            6,
            6.25,
            6.5,
            6.75,
            7,
            7.25,
            7.5,
            7.75,
            8,
            8.25,
            8.5,
            8.75,
            9,
            9.25,
            9.5,
            9.75,
            10,
            10.5,
            10.75,
            11.5,
            11.75,
        ]
        CL = [
            -0.0105,
            0.0047,
            0.0103,
            0.0386,
            0.0702,
            0.0884,
            0.1103,
            0.1289,
            0.1542,
            0.1897,
            0.224,
            0.2714,
            0.3355,
            0.4065,
            0.4772,
            0.8569,
            0.9316,
            0.9834,
            1.0275,
            1.0711,
            1.1236,
            1.1483,
            1.1729,
            1.1995,
            1.2338,
            1.2672,
            1.2977,
            1.3313,
            1.365,
            1.3918,
            1.4219,
            1.4531,
            1.4823,
            1.5037,
            1.5223,
            1.5379,
            1.5728,
            1.6028,
            1.6248,
            1.649,
            1.6752,
            1.7032,
            1.7335,
            1.7484,
            1.7673,
            1.7899,
            1.8152,
            1.8443,
            1.8582,
            1.8685,
            1.886,
            1.911,
            1.9414,
            1.9463,
            1.9457,
            1.9623,
            1.9937,
            2.022,
            1.9896,
            1.9883,
            2.0272,
            2.0814,
            1.9821,
            1.9679,
            2.0778,
            1.8149,
            1.4189,
        ]

        CLalpha = InterpolatedUnivariateSpline(np.radians(alpha_data), CL)
        a0 = CLalpha.derivative(1)

        alpha0 = -np.radians(5)
        A = np.pi * R**2

        r2R = np.arange(1, n_elements + 1) / n_elements  # Local radius to total rotor

        c2R = c / R

        M_local = (r2R) * (omega * R / V_m)  # Local Mach number (omega*R = V_tip which is constant)

        a = a0(takeoffAOA) / (
            np.sqrt(1 - M_local**2)
        )  # Lift curve slope corrected for mach number

        def v12Omegar(theta):
            return (
                (a * b * c2R)
                / (16 * np.pi * r2R)
                * (-1 + np.sqrt(1 + (32 * np.pi * theta * r2R) / (a * b * c2R)))
            )

        alpha = takeoffAOA * np.ones(n_elements)

        def func(x):
            return alpha - (x - np.arctan(v12Omegar(x)))

        theta = scipy.optimize.fsolve(
            func,
            0.1
            * np.ones(
                n_elements,
            ),
        )

        # plt.plot(r2R, np.degrees(theta))
        # plt.title("Required twist angle")
        # plt.show()

        # S1223
        cl = CLalpha(alpha)
        cd = cl / 50
        DctDr2R = b * r2R**2 * c2R * cl / (2 * np.pi)
        # Create spline of data points in order
        funCT = InterpolatedUnivariateSpline(r2R, DctDr2R)
        Ct_lossless = funCT.integral(x0, 1)
        # Loss of lift from the tip of the rotor
        if Ct_lossless < 0.006:
            B = 1 - 0.06 / b
        else:
            B = 1 - np.sqrt(2.27 * Ct_lossless - 0.01) / b

        Ct = Ct_lossless - funCT.integral(B, 1)

        Dcq0Dr2R = b * r2R**3 * c2R * cd / (2 * np.pi)
        funCQ0 = InterpolatedUnivariateSpline(r2R, Dcq0Dr2R)
        CQ_profile = funCQ0.integral(x0, 1)
        DcqiDr2R = b * r2R**3 * c2R * cl * v12Omegar(theta) / (2 * np.pi)
        funCQi = InterpolatedUnivariateSpline(r2R, DcqiDr2R)
        CQ_induced = funCQi.integral(B, 1)

        DCQ_I = 0.01 * CQ_induced

        Cq = (CQ_profile + CQ_induced + DCQ_I) / 0.95
        torque = rho * A * (omega * R) ** 2 * Cq

        T = rho * A * (omega * R) ** 2 * Ct
        if coaxial:
            T *= 0.88
        pow = rho * A * V_tip**3 * Cq / 550 / 1.341 * 1000
        pow = power(T, R)
    sigma = b * c / (np.pi * R)
    # print(Ct/sigma)
    print(torque)

    # Rotor Weight:
    x_cord_top = np.flip(
        [
            1,
            0.99838,
            0.99417,
            0.98825,
            0.98075,
            0.97111,
            0.95884,
            0.94389,
            0.92639,
            0.90641,
            0.88406,
            0.85947,
            0.83277,
            0.80412,
            0.77369,
            0.74166,
            0.70823,
            0.6736,
            0.63798,
            0.60158,
            0.56465,
            0.52744,
            0.49025,
            0.4534,
            0.41721,
            0.38193,
            0.34777,
            0.31488,
            0.28347,
            0.2537,
            0.22541,
            0.19846,
            0.17286,
            0.14863,
            0.12591,
            0.10482,
            0.08545,
            0.06789,
            0.05223,
            0.03855,
            0.02694,
            0.01755,
            0.01028,
            0.00495,
            0.00155,
            0.00005,
        ]
    )
    x_chord_bottom = [
        0.00005,
        0.00044,
        0.00264,
        0.00789,
        0.01718,
        0.03006,
        0.04627,
        0.06561,
        0.08787,
        0.11282,
        0.1402,
        0.17006,
        0.20278,
        0.2384,
        0.27673,
        0.3175,
        0.36044,
        0.40519,
        0.45139,
        0.4986,
        0.54639,
        0.59428,
        0.64176,
        0.68832,
        0.73344,
        0.7766,
        0.81729,
        0.855,
        0.88928,
        0.91966,
        0.94573,
        0.96693,
        0.98255,
        0.99268,
        0.99825,
        1,
    ]
    y_cord_top = np.flip(
        [
            0,
            0.00126,
            0.00494,
            0.01037,
            0.01646,
            0.0225,
            0.02853,
            0.03476,
            0.04116,
            0.04768,
            0.05427,
            0.06089,
            0.06749,
            0.07402,
            0.08044,
            0.08671,
            0.09277,
            0.09859,
            0.10412,
            0.10935,
            0.11425,
            0.11881,
            0.12303,
            0.12683,
            0.13011,
            0.13271,
            0.13447,
            0.13526,
            0.13505,
            0.13346,
            0.13037,
            0.12594,
            0.12026,
            0.11355,
            0.10598,
            0.0977,
            0.08879,
            0.0794,
            0.06965,
            0.05968,
            0.04966,
            0.03961,
            0.02954,
            0.01969,
            0.01033,
            0.00178,
        ]
    )
    y_chord_bottom = [
        0.00178,
        -0.00561,
        -0.0112,
        -0.01427,
        -0.0155,
        -0.01584,
        -0.01532,
        -0.01404,
        -0.01202,
        -0.00925,
        -0.00563,
        -0.00075,
        0.00535,
        0.01213,
        0.01928,
        0.02652,
        0.03358,
        0.04021,
        0.04618,
        0.05129,
        0.05534,
        0.0582,
        0.05976,
        0.05994,
        0.05872,
        0.05612,
        0.05219,
        0.04706,
        0.04088,
        0.03387,
        0.02624,
        0.01822,
        0.0106,
        0.00468,
        0.00115,
        0,
    ]

    S1223_top = InterpolatedUnivariateSpline(x_cord_top, y_cord_top)
    S1223_bottom = InterpolatedUnivariateSpline(x_chord_bottom, y_chord_bottom)

    Area_top = S1223_top.integral(0, 1)
    Area_bot = S1223_bottom.integral(0, 1)
    Area = (Area_top - Area_bot) * c

    Rotor_mass = b * const["bladeDensity"] * R * Area * const["fillFactor"]

    if print_results:
        print("Given the mass of " + str(M) + "kg, the following applied: \n")
        print(f"Rotor Thrust: {T}[N]")
        print("Rotor Radius: " + str(R) + "[m]")
        print("Rotor Power: " + str(pow) + "[W]")
        print("Total Aircraft power: " + str(pow * N_rotors) + "[W]")
        print("Weight per rotor with " + str(b) + " blades is " + str(Rotor_mass) + "[kg]")
        print(f"Total Rotor mass for entire aircraft: {Rotor_mass*N_rotors}[kg]")

    return R, T, pow / 745.7, pow * N_rotors, Rotor_mass * N_rotors


def RotorInCruise(V_cr, omega, R):
    r = np.linspace(0, R, 100)
    V_total = np.sqrt((V_cr) ** 2 + omega**2 * r**2)
    phi = np.arctan(V_cr / (omega * r))

    # Forward component of lift:

    plt.plot(r, np.cos(phi))
    plt.show()


if __name__ == "__main__":
    gm = const["gravityMars"]
    rho = const["airDensity"]  # kg/m3
    MTOW = 3000 * gm

    # Blade Element Method with ideal twist
    N_rotor = 4

    RadiusMassElementMomentum(3000, 4, 6, coaxial=True, V_tip=168)

    cruiseSpeed = np.linspace(112, 154)
    a, b = max_tipSpeed(cruiseSpeed)
    V_tip = a * cruiseSpeed + b

    R = np.empty(np.shape(V_tip))
    for i in range(len(V_tip)):
        R[i] = RadiusMassElementMomentum(3000, 4, 6, coaxial=True, V_tip=V_tip[i])[0]

    def oneOverX(x, a, b, c, d):
        return a / (b * x + c) + d

    popt2, pcov2, infoDict, msg, ier = curve_fit(oneOverX, V_tip, R, full_output=True)
    print(f"Rotor Radius =  {popt2[0]}/({popt2[1]}x + {popt2[2]}) + {popt2[3]}")
    print(pcov2)
    plt.scatter(V_tip, R, alpha=0.5)
    plt.plot(V_tip, oneOverX(V_tip, *popt2), "r--")
    plt.xlabel("Tip speed [m/s]")
    plt.ylabel("Rotor radius [m]")
    plt.show()

    plt.plot(
        cruiseSpeed,
        oneOverX(a * cruiseSpeed + b, *popt2),
        label=f"{np.round(popt2[0], 2)}/"
        f"({np.round(popt2[1]*a, 2)}cruiseSpeed + "
        f"{np.round(popt2[1]*b + popt2[2], 2)}) + "
        f"{np.round(popt2[3], 2)})",
    )
    plt.xlabel("Cruise speed [m/s]")
    plt.ylabel("Rotor radius [m]")
    plt.legend()
    plt.show()
    print(
        f"rotorRadius = {np.round(popt2[0], 2)}/({np.round(popt2[1]*a, 2)}cruiseSpeed + "
        f"{np.round(popt2[1]*b + popt2[2], 2)}) + {np.round(popt2[3], 2)})"
    )

    """'
    hp = T*v1*0.01315
    DL = T/A
    
    plt.scatter(R, T)
    plt.plot(R, T_init*np.ones(np.shape(R)))
    plt.show()
    
    #Using above values:
    
    R1=21.5#[m]
    A1=np.pi* R1**2
    omega1 = 200/R1
    c1=R1/25
    sigma = 6*c1/(np.pi*R1)
    Ct = T_init/(rho*A1*(omega1*R1)**2)
    
    T=np.empty(23,)
    hp=np.empty(23,)
    DL=np.empty(23,)
    for R in range(2, 25):
        c = R / 25
        A = np.pi * R ** 2
        omega = 200 / R
        v1 = np.sqrt(T_init / (2 * rho * A))
        theta_tip = 10 * np.pi / 180
        T[R-2] = b * rho * omega ** 2 * a * c / 2 * (
                    theta_tip * R ** 3 / 2 - (v1 * (np.log(v1 ** 2 + (omega * R) ** 2) - np.log(v1 ** 2)) / (2 * omega)))
        hp[R-2] = T[R-2] * v1 * 0.01315
        DL[R-2] = T[R-2] / A
    
    R=np.arange(2, 25, 1)
    
    plt.scatter(R, T)
    plt.plot(R, T_init*np.ones(np.shape(R)))
    plt.show()
    
    plt.scatter(R, hp)
    plt.show()
    
    plt.scatter(R, DL)
    plt.show()
    """
