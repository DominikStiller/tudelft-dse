import math
import warnings

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from tqdm import tqdm

warnings.filterwarnings("ignore")  # TODO: DELETE THIS!!!

const = {
    "gravityMars": 3.71,
    "airDensity": 0.01,
    "soundSpeed": 220,
    "cl": 1.6,
    "cd": 0.03,
    "irradiance": 590,
    "solarPanelDensity": 1.76,
    "takeoffBatteryPowerDensity": 1317,
    "takeoffBatteryEnergyDensity": 437,
    "bladeDensity": 1500,
    "allowedStress": 1.2 * 10**9,
    "t/c": 0.012,
    "visc_cr": 5.167e-4,
    "cruiseSpeed": 400 / 3.6,
    "batteryVolume": 450,  # Wh/L
    "ultimateLoad": 1.5,
    "takeOffLoad": 1.1,
    "designRange": 1e6,
    "takeOffTime": 300,
    "margin": 0.0,
    "payloadMass": 350,
    "maxMass": 2700,
    "oswald": 0.9,
    "fillFactor": 0.084,
}

rotorParameters = {
    "maxEnginePower": 52000,
    "Radius": 8.2,
    "chord": 8.2 / 30,
    "cutout": 0.15,
    "N_blades": 8,
    "coaxial": True,
    "CL_data": [
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
    ],
    "CD_data": [
        0.10503,
        0.10263,
        0.10225,
        0.09558,
        0.09134,
        0.08839,
        0.08464,
        0.08072,
        0.07715,
        0.0726,
        0.06862,
        0.06377,
        0.05889,
        0.05411,
        0.04983,
        0.02451,
        0.02289,
        0.02207,
        0.02169,
        0.02145,
        0.02167,
        0.02246,
        0.02337,
        0.02372,
        0.02405,
        0.02443,
        0.02476,
        0.02513,
        0.02561,
        0.02601,
        0.02644,
        0.0269,
        0.02744,
        0.02784,
        0.02811,
        0.02823,
        0.02903,
        0.02996,
        0.03079,
        0.03161,
        0.0324,
        0.03323,
        0.03427,
        0.03535,
        0.03642,
        0.03739,
        0.03826,
        0.03914,
        0.04061,
        0.04213,
        0.04332,
        0.04411,
        0.04485,
        0.04684,
        0.04896,
        0.05011,
        0.05041,
        0.0513,
        0.05513,
        0.05721,
        0.05688,
        0.05612,
        0.06371,
        0.06835,
        0.06333,
        0.09405,
        0.16839,
    ],
    "alpha_data": np.radians(
        [
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
    ),
    "takeoff_AoA": np.radians(8),
    "cruise_AoA": np.radians(3.5),
    "alpha_0": np.radians(-5),
    "blade_twist": [
        14.32000402,
        10.29123294,
        8.32003476,
        7.12117696,
        6.31005798,
        5.72821286,
        5.29796352,
        4.97731679,
        4.74279731,
        4.58226976,
        4.49233726,
        4.47898847,
        4.56296873,
        4.79802654,
        5.34554722,
    ],
    "Rear_blade_twist": 0,
    "n_elements": 20,
}


def TakeoffRotor(V_tip, prnt, n_elements=rotorParameters["n_elements"]):
    # Rename constants
    gm = const["gravityMars"]
    rho = const["airDensity"]
    V_m = const["soundSpeed"]
    R = rotorParameters["Radius"]
    b = rotorParameters["N_blades"]
    v_tip = V_tip
    coaxial = rotorParameters["coaxial"]
    c = rotorParameters["chord"]
    x0 = rotorParameters["cutout"]
    omega = v_tip / R

    takeoffAOA = rotorParameters["takeoff_AoA"]

    CLalpha = InterpolatedUnivariateSpline(
        rotorParameters["alpha_data"], rotorParameters["CL_data"]
    )
    CDalpha = InterpolatedUnivariateSpline(
        rotorParameters["alpha_data"], rotorParameters["CD_data"]
    )
    a0 = CLalpha.derivative(1)

    A = np.pi * R**2

    r2R = np.arange(1, n_elements + 1) / n_elements  # Local radius to total rotor

    c2R = c / R / r2R
    c2R[c2R * R > 0.75] = 0.75 / R

    M_local = (r2R) * (omega * R / V_m)  # Local Mach number (omega*R = V_tip which is constant)

    a = a0(takeoffAOA) / (np.sqrt(1 - M_local**2))  # Lift curve slope corrected for mach number
    V_c = 0

    def v12Omegar(theta):
        v_1cr = (
            -(omega / 2 * a * c2R * R * b + 4 * np.pi * V_c)
            + np.sqrt(
                (omega / 2 * a * c2R * R * b + 4 * np.pi * V_c) ** 2
                + 8
                * np.pi
                * b
                * omega**2
                * a
                * c2R
                * R
                * r2R
                * R
                * (theta - V_c / (omega * r2R * R))
            )
        ) / (8 * np.pi * omega * r2R * R)
        v_1 = ((a * b * c2R) / (16 * np.pi * r2R)) * (
            -1 + np.sqrt(1 + (32 * np.pi * theta * r2R) / (a * b * c2R))
        )

        return v_1

    alpha = takeoffAOA * np.ones(n_elements)

    def func(x):
        return alpha - (x - np.arctan(v12Omegar(x)))

    theta = scipy.optimize.fsolve(
        func,
        0.5
        * np.ones(
            n_elements,
        ),
    )
    rotorParameters["blade_twist"] = theta + rotorParameters["alpha_0"]
    # plt.plot(r2R, func(theta))
    # plt.show()

    # S1223
    cl = CLalpha(func(theta) + alpha)
    cd = CDalpha(func(theta) + alpha)
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
    CQ_induced = funCQi.integral(x0, B)

    Pratiofunc = r2R**3 * (1 - np.sqrt(1 - (2 * Ct / (r2R**2)))) ** 2
    DCQ_I = 0
    if any(np.isfinite(Pratiofunc)) == True:
        ii = np.isfinite(Pratiofunc)
        funcPratio = InterpolatedUnivariateSpline(r2R[ii], Pratiofunc[ii])

        DCQ_I = 1 / Ct * funcPratio.integral(np.sqrt(2 * Ct), 1) * CQ_induced

    Cq = (CQ_profile + CQ_induced + DCQ_I) / 0.95
    Vblade = omega * R
    T = rho * A * (Vblade) ** 2 * Ct

    pow = rho * A * V_tip**3 * Cq
    torque = rho * A * (omega * R) ** 2 * Cq
    # pow = power(T, R)
    sigma = b * c / (np.pi * R)
    ct2sigma = Ct / sigma
    cq2sigma = Cq / sigma

    Vinf = 0
    PitchDiam = r2R * np.pi * np.arctan(rotorParameters["blade_twist"])
    Ftilde = T * 2.64 / (R - R * x0)
    rHAT = (r2R * R - x0 * R) / (R - x0 * R)
    rHAT = rHAT.clip(min=0)
    m = 1
    a_f = 1
    n = 0.2
    fx = Ftilde * rHAT**m * (1 - rHAT / a_f) ** n
    ftheta = fx / np.pi * PitchDiam / r2R
    Vix = (
        np.sqrt(Vinf**2 / 4 + fx / (4 * rho * np.pi * r2R * R)) - Vinf / 2
    )  # Axial Flow Velocity
    Vitheta = ftheta / (2 * rho * np.pi * (Vix + Vinf))  # Tangential Flow Velocity
    Vitheta = np.nan_to_num(Vitheta, nan=0)

    AddHP, collective = AdditionalClimbPower(T, 2)
    v = const["visc_cr"]
    Re_vortex = 2 * omega * R * c2R[-1] * R / v * Ct / sigma
    if prnt == True:
        print(f"###########UPPER ROTOR############")
        print(f"Using a rotor radius of: {R}[m]")
        print(f"Cutout: {x0}")
        print(f"Chord: {R*c2R}")
        print(f"Re at tip: {Re_vortex}")
        print(f"Downwash: {np.max(Vix)}[m/s] (AXIAL) and {np.max(Vitheta)}[m/s](Tangential)")
        print(f"Front Blade Twist: {np.degrees(rotorParameters['blade_twist'])}")
        print(f"Hover Thrust: {T}[N]")
        print(f"Hover Torque: {torque}[Nm]")
        print(f"Hover Power: {pow / 1000}[kW]")
        print(f"Vertical Climb Power per rotor: {(pow+AddHP)/1000}[kW]")
        print(f"Change in collective for climb: {np.degrees(collective)}")
        print(f"Hover RPM: {V_tip/R * 60/(2*np.pi)}")
    if any(alpha > np.radians(11.5)):
        T = 0
        torque = 0
        pow = 0
        AddHP = 0
    return T, torque, pow, AddHP


def SecondRotor(
    vtipmax, thetaRotor1, T_Rotor1, n_elements=rotorParameters["n_elements"], prnt=False
):
    gm = const["gravityMars"]
    rho = const["airDensity"]
    V_m = const["soundSpeed"]
    R = rotorParameters["Radius"]
    b = rotorParameters["N_blades"]
    coaxial = rotorParameters["coaxial"]
    C = rotorParameters["chord"]
    x0 = rotorParameters["cutout"]

    takeoffAOA = rotorParameters["cruise_AoA"]

    CLalpha = InterpolatedUnivariateSpline(
        rotorParameters["alpha_data"], rotorParameters["CL_data"]
    )
    CDalpha = InterpolatedUnivariateSpline(
        rotorParameters["alpha_data"], rotorParameters["CD_data"]
    )
    a0 = CLalpha.derivative(1)

    A = np.pi * R**2

    r2R = np.arange(1, n_elements + 1) / n_elements  # Local radius to total rotor

    c2R = C / R / r2R
    c2R[c2R * R > 0.75] = 0.75 / R

    # Induced Swirl from previous rotor
    Vinf = 0
    PitchDiam = r2R * np.pi * np.arctan(rotorParameters["blade_twist"])
    Ftilde = T_Rotor1 * 2.64 / (R - R * x0)
    rHAT = (r2R * R - x0 * R) / (R - x0 * R)
    rHAT = rHAT.clip(min=0)
    m = 1
    a_f = 1
    n = 0.2
    fx = Ftilde * rHAT**m * (1 - rHAT / a_f) ** n
    ftheta = fx / np.pi * PitchDiam / r2R
    Vix = (
        np.sqrt(Vinf**2 / 4 + fx / (4 * rho * np.pi * r2R * R)) - Vinf / 2
    )  # Axial Flow Velocity
    Vitheta = ftheta / (2 * rho * np.pi * (Vix + Vinf))  # Tangential Flow Velocity
    Vitheta = np.nan_to_num(Vitheta, nan=0)
    Vitheta[-1] = Vitheta[-2]

    def funcomega(ome):
        return vtipmax - np.sqrt(np.average(Vix[Vix > 0]) ** 2 + (ome * R + Vitheta[-2]) ** 2)

    omega = scipy.optimize.fsolve(funcomega, np.array([15]))[0]
    omega_actual = omega

    M_local = (r2R) * (vtipmax / V_m)  # Local Mach number (omega*R = V_tip which is constant)

    a = a0(takeoffAOA) / (np.sqrt(1 - M_local**2))  # Lift curve slope corrected for mach number

    alpha = takeoffAOA * np.ones(n_elements)
    V_c = np.average(Vix[Vix > 0])

    def v12Omegar(x):
        v_1cr = (
            -(omega_actual / 2 * a * R * c2R * b + 4 * np.pi * V_c)
            + np.sqrt(
                (omega_actual / 2 * a * R * c2R * b + 4 * np.pi * V_c) ** 2
                + 8
                * np.pi
                * b
                * (omega_actual**2)
                * a
                * R
                * c2R
                * R
                * r2R
                * (x - V_c / (omega_actual * R * r2R))
            )
        ) / (8 * np.pi)

        return v_1cr

    def func(x):
        return x - (alpha + np.arctan((v12Omegar(x) + V_c) / (omega * R * r2R + Vitheta)))

    i = 0
    theta = [0, 0]
    while theta[0] == theta[-1]:
        theta = scipy.optimize.fsolve(
            func,
            i
            * np.ones(
                n_elements,
            ),
        )
        i += 0.005

    rotorParameters["Rear_blade_twist"] = theta + rotorParameters["alpha_0"]

    # S1223
    cl = CLalpha(func(theta) + alpha)
    cd = CDalpha(func(theta) + alpha)
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
    DcqiDr2R = (
        b
        * r2R**3
        * c2R
        * cl
        * ((v12Omegar(theta) + V_c) / (omega * R * r2R + Vitheta))
        / (2 * np.pi)
    )
    funCQi = InterpolatedUnivariateSpline(r2R, DcqiDr2R)
    CQ_induced = funCQi.integral(x0, B)

    Pratiofunc = r2R**3 * (1 - np.sqrt(1 - (2 * Ct / (r2R**2)))) ** 2
    DCQ_I = 0
    if any(np.isfinite(Pratiofunc)) == True:
        ii = np.isfinite(Pratiofunc)
        funcPratio = InterpolatedUnivariateSpline(r2R[ii], Pratiofunc[ii])

        DCQ_I = 1 / Ct * funcPratio.integral(np.sqrt(2 * Ct), 1) * CQ_induced

    Cq = (CQ_profile + CQ_induced + DCQ_I) / 0.95

    T = rho * A * (vtipmax) ** 2 * Ct
    pow = rho * A * (vtipmax) ** 3 * Cq
    torque = rho * A * (vtipmax) ** 2 * Cq
    sigma = b * C / (np.pi * R)
    ct2sigma = Ct / sigma
    cq2sigma = Cq / sigma

    AddHP, collective = AdditionalClimbPower(T, V_c)
    Vinf2 = Vix
    PitchDiam2 = r2R * np.pi * np.arctan(rotorParameters["Rear_blade_twist"])
    Ftilde2 = T * 2.64 / (R - R * x0)
    rHAT2 = (r2R * R - x0 * R) / (R - x0 * R)
    rHAT2 = rHAT2.clip(min=0)
    m = 1
    a_f = 1
    n = 0.2
    fx2 = Ftilde2 * rHAT2**m * (1 - rHAT2 / a_f) ** n
    ftheta2 = fx2 / np.pi * PitchDiam2 / r2R
    Vix2 = (
        np.sqrt(Vinf2**2 / 4 + fx2 / (4 * rho * np.pi * r2R * R)) - Vinf2 / 2
    )  # Axial Flow Velocity
    Vitheta2 = ftheta2 / (2 * rho * np.pi * (Vix2 + Vinf2))  # Tangential Flow Velocity
    Vitheta2 = np.nan_to_num(Vitheta2, nan=0)

    v = const["visc_cr"]
    Re_vortex = 2 * omega * R * c2R[-1] * R / v * Ct / sigma

    if prnt == True:
        print(f"##############REAR BLADE##############")
        print(f"Using a rotor radius of: {R}[m]")
        print(f"Cutout: {x0}")
        print(f"Chord: {R * c2R}")
        print(f"Re at tip: {Re_vortex}")
        print(f"Rear Blade Twist: {np.degrees(rotorParameters['Rear_blade_twist'])}")
        print(f"Hover Thrust: {T}[N]")
        print(f"Hover Torque: {torque}[Nm]")
        print(f"Hover Power: {pow / 1000}[kW]")
        print(f"Hover RPM: {omega_actual * 60/(2*np.pi)}")
        print(
            f"Downwash: {np.max(Vix2+Vix)}[m/s] (AXIAL) and {np.max(abs(Vitheta2))}[m/s](Tangential)"
        )
        print(f"Vertical Climb Power: {(pow+AddHP)/1000}[kW]")
        print(f"Change in collective for climb: {np.degrees(collective)}")

    return T, torque, pow, AddHP


def DetermineCutout(MTOM, tipspeed, plot=False):
    T = 10000 * 3.71 * np.ones(100)
    Torque = np.zeros(100)
    i = 0
    c, d = 0, 0
    while not min(T) <= (MTOM * 3.71):
        rotorParameters["cutout"] += 0.01
        T_1, Torque_1, a, b = TakeoffRotor(
            tipspeed, prnt=False, n_elements=rotorParameters["n_elements"]
        )
        T_2, Torque_2, c, d = SecondRotor(
            tipspeed,
            rotorParameters["blade_twist"],
            T_1,
            n_elements=rotorParameters["n_elements"],
            prnt=False,
        )

        T[i] = 2 * (T_1 + T_2)
        Torque[i] = 2 * (Torque_1 + Torque_2)
        i += 1

    def round_to_nearest_multiple(number):
        multiple = 1 / 20
        rounded_number = round(number / multiple) * multiple
        return rounded_number

    rotorParameters["cutout"] = math.floor(
        rotorParameters["cutout"] / (1 / rotorParameters["n_elements"])
    ) * (1 / rotorParameters["n_elements"])

    T = T[:i]
    Torque = Torque[:i]
    if plot == True:
        plt.scatter(np.linspace(0, rotorParameters["cutout"], len(T)), T)
        plt.show()
        plt.scatter(np.linspace(0, rotorParameters["cutout"], len(T)), Torque)
        plt.show()
    print(f"Max possible cutout to lift a mass of {MTOM}[kg] is {rotorParameters['cutout']*100}%")


def OptimizeCruise(Drag, Vcr, plot=False):
    def FirstRotorThrustAndPower(x, V_c, aoa, plot=False):
        gm = const["gravityMars"]
        rho = const["airDensity"]
        V_m = const["soundSpeed"]
        R = rotorParameters["Radius"]
        b = rotorParameters["N_blades"]
        coaxial = rotorParameters["coaxial"]
        c = rotorParameters["chord"]
        x0 = rotorParameters["cutout"]
        n_elements = rotorParameters["n_elements"]

        takeoffAOA = aoa

        theta0 = rotorParameters["blade_twist"]
        omega, collective = x

        CLalpha = InterpolatedUnivariateSpline(
            rotorParameters["alpha_data"], rotorParameters["CL_data"]
        )
        CDalpha = InterpolatedUnivariateSpline(
            rotorParameters["alpha_data"], rotorParameters["CD_data"]
        )
        a0 = CLalpha.derivative(1)

        A = np.pi * R**2

        r2R = np.arange(1, n_elements + 1) / n_elements  # Local radius to total rotor

        c2R = c / R
        c2R = c / R / r2R
        c2R[c2R * R > 0.75] = 0.75 / R

        theta = theta0 + collective - rotorParameters["alpha_0"]

        # def funcomega(ome):
        # return vtipmax - np.sqrt(V_c ** 2 + (ome * R) ** 2)

        # omega = scipy.optimize.fsolve(funcomega, np.array([2]))[0]
        omega_actual = omega
        #
        # M_local = (r2R) * (vtipmax / V_m)  # Local Mach number (omega*R = V_tip which is constant)

        # a = a0(takeoffAOA) / (np.sqrt(1 - M_local ** 2))  # Lift curve slope corrected for mach number
        a = a0(takeoffAOA)

        def v12Omegar(x):
            v_1cr = (
                -(omega_actual / 2 * a * R * c2R * b + 4 * np.pi * V_c)
                + np.sqrt(
                    (omega_actual / 2 * a * R * c2R * b + 4 * np.pi * V_c) ** 2
                    + 8
                    * np.pi
                    * b
                    * (omega_actual**2)
                    * a
                    * R
                    * c2R
                    * R
                    * r2R
                    * (x - V_c / (omega_actual * R * r2R))
                )
            ) / (8 * np.pi)

            return v_1cr

        vtipmax = np.sqrt((V_c + v12Omegar(theta)[-1]) ** 2 + (omega * R))
        alpha = theta - np.arctan((v12Omegar(theta) + V_c) / (omega * r2R * R))
        # S1223
        cl = CLalpha(alpha)
        cd = CDalpha(alpha)
        if any(alpha > np.radians(11)):
            cl[alpha > np.radians(11)] = 0.5 * np.cos(alpha[alpha > np.radians(11)])
            cd[alpha > np.radians(11)] = 1.28 * np.sin(alpha[alpha > np.radians(11)])
        if any(alpha < np.radians(-5)):
            cl[alpha < np.radians(-5)] = -0.000001
            cd[alpha < np.radians(-5)] = 1.28 * np.sin(alpha[alpha < np.radians(-5)])

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

        DcqiDr2R = (
            b * r2R**3 * c2R * cl * (v12Omegar(theta) + V_c) / (omega * r2R * R) / (2 * np.pi)
        )

        funCQi = InterpolatedUnivariateSpline(r2R, DcqiDr2R)
        CQ_induced = funCQi.integral(x0, B)

        Pratiofunc = r2R**3 * (1 - np.sqrt(1 - (2 * Ct / (r2R**2)))) ** 2
        DCQ_I = 0
        if any(np.isfinite(Pratiofunc)) == True:
            ii = np.isfinite(Pratiofunc)
            funcPratio = InterpolatedUnivariateSpline(r2R[ii], Pratiofunc[ii])

            DCQ_I = 1 / Ct * funcPratio.integral(np.sqrt(2 * Ct), 1) * CQ_induced

        Cq = (CQ_profile + CQ_induced + DCQ_I) / 0.95

        T = rho * A * vtipmax**2 * Ct
        Torque = rho * A * vtipmax**2 * Cq
        pow = rho * A * vtipmax**3 * Cq
        sigma = b * c / (np.pi * R)
        pow += AdditionalClimbPower(T, const["cruiseSpeed"])[0]

        if plot == True:
            print(f"Rotor 1 AOA: {np.degrees(alpha)}")
            plt.plot(r2R, funCQi(r2R))
            plt.plot(r2R, funCQ0(r2R))
            plt.plot(r2R, funcPratio(r2R))
            plt.show()
        # if any(alpha > np.radians(11.5)):
        #     T = 0
        #     torque = 0
        #     pow = 0
        #     AddHP = 0
        return T, pow, Torque

    def SecondRotorThrustAndPower(x, T_1, V_c, aoa, plot=False):
        gm = const["gravityMars"]
        rho = const["airDensity"]
        V_m = const["soundSpeed"]
        R = rotorParameters["Radius"]
        b = rotorParameters["N_blades"]
        coaxial = rotorParameters["coaxial"]
        c = rotorParameters["chord"]
        x0 = rotorParameters["cutout"]
        n_elements = rotorParameters["n_elements"]
        takeoffAOA = aoa
        theta0 = rotorParameters["Rear_blade_twist"]
        CLalpha = InterpolatedUnivariateSpline(
            rotorParameters["alpha_data"], rotorParameters["CL_data"]
        )
        CDalpha = InterpolatedUnivariateSpline(
            rotorParameters["alpha_data"], rotorParameters["CD_data"]
        )
        a0 = CLalpha.derivative(1)

        A = np.pi * R**2

        r2R = np.arange(1, n_elements + 1) / n_elements  # Local radius to total rotor

        c2R = c / R
        c2R = c / R / r2R
        c2R[c2R * R > 0.75] = 0.75 / R

        omega, collective = x

        # Induced Swirl from previous rotor
        Vinf = V_c
        PitchDiam = r2R * np.pi * np.arctan(rotorParameters["blade_twist"])
        Ftilde = T_1 * 2.64 / (R - R * x0)
        rHAT = (r2R * R - x0 * R) / (R - x0 * R)
        rHAT = rHAT.clip(min=0)
        m = 1
        a_f = 1
        n = 0.2
        fx = Ftilde * rHAT**m * (1 - rHAT / a_f) ** n
        ftheta = fx / np.pi * PitchDiam / r2R
        Vix = (
            np.sqrt(Vinf**2 / 4 + fx / (4 * rho * np.pi * r2R * R)) - Vinf / 2
        )  # Axial Flow Velocity
        Vitheta = ftheta / (2 * rho * np.pi * (Vix + Vinf))  # Tangential Flow Velocity
        Vitheta = np.nan_to_num(Vitheta, nan=0)
        Vitheta[-1] = Vitheta[-2]

        # def funcomega(ome):
        #     return vtipmax - np.sqrt(Vix[-2] ** 2 + (ome * R + Vitheta[-2]) ** 2)
        #
        # omega = scipy.optimize.fsolve(funcomega, np.array([2]))[0]
        omega_actual = omega

        # M_local = (r2R) * (vtipmax / V_m)  # Local Mach number (omega*R = V_tip which is constant)
        #
        # a = a0(takeoffAOA) / (np.sqrt(1 - M_local ** 2))  # Lift curve slope corrected for mach number
        a = a0(takeoffAOA)

        alpha = takeoffAOA * np.ones(n_elements)
        V_c += np.average(Vix[Vix > 0])

        def v12Omegar(x):
            v_1cr = (
                -(omega_actual / 2 * a * R * c2R * b + 4 * np.pi * V_c)
                + np.sqrt(
                    (omega_actual / 2 * a * R * c2R * b + 4 * np.pi * V_c) ** 2
                    + 8
                    * np.pi
                    * b
                    * (omega_actual**2)
                    * a
                    * R
                    * c2R
                    * R
                    * r2R
                    * (x - V_c / (omega_actual * R * r2R))
                )
            ) / (8 * np.pi)

            return v_1cr

        theta = theta0 + collective - rotorParameters["alpha_0"]
        alpha = theta - np.arctan((v12Omegar(theta) + V_c) / (omega * r2R * R + Vitheta))
        vtipmax = np.sqrt(
            (omega * R + Vitheta[-2]) ** 2
            + (V_c + np.average(Vix[Vix > 0]) + v12Omegar(theta)[-1]) ** 2
        )
        # S1223
        cl = CLalpha(alpha)
        cd = CDalpha(alpha)
        if any(alpha > np.radians(11)):
            cl[alpha > np.radians(11)] = 0.5 * np.cos(alpha[alpha > np.radians(11)])
            cd[alpha > np.radians(11)] = 1.28 * np.sin(alpha[alpha > np.radians(11)])
        if any(alpha < np.radians(-5)):
            cl[alpha < np.radians(-5)] = -0.000001
            cd[alpha < np.radians(-5)] = 1.28 * np.sin(alpha[alpha < np.radians(-5)])
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

        DcqiDr2R = (
            b
            * r2R**3
            * c2R
            * cl
            * (v12Omegar(theta) + V_c)
            / (omega * r2R * R + Vitheta)
            / (2 * np.pi)
        )

        funCQi = InterpolatedUnivariateSpline(r2R, DcqiDr2R)
        CQ_induced = funCQi.integral(x0, B)

        Pratiofunc = r2R**3 * (1 - np.sqrt(1 - (2 * Ct / (r2R**2)))) ** 2
        DCQ_I = 0
        funcPratio = r2R * 0
        if any(np.isfinite(Pratiofunc)) == True:
            ii = np.isfinite(Pratiofunc)
            funcPratio = InterpolatedUnivariateSpline(r2R[ii], Pratiofunc[ii])

            DCQ_I = 1 / Ct * funcPratio.integral(np.sqrt(2 * Ct), 1) * CQ_induced

        Cq = (CQ_profile + CQ_induced + DCQ_I) / 0.95

        T = rho * A * vtipmax**2 * Ct
        Torque = rho * A * vtipmax**2 * Cq
        pow = rho * A * vtipmax**3 * Cq
        sigma = b * c / (np.pi * R)
        pow += AdditionalClimbPower(T, const["cruiseSpeed"])[0]

        Vinf2 = Vcr + Vix
        PitchDiam2 = r2R * np.pi * np.arctan(theta)
        Ftilde2 = T_1 * 2.64 / (R - R * x0)
        rHAT2 = (r2R * R - x0 * R) / (R - x0 * R)
        rHAT2 = rHAT2.clip(min=0)
        m = 1
        a_f = 1
        n = 0.2
        fx2 = Ftilde2 * rHAT2**m * (1 - rHAT2 / a_f) ** n
        ftheta2 = fx2 / np.pi * PitchDiam2 / r2R
        Vix2 = np.sqrt(Vinf2**2 / 4 + fx / (4 * rho * np.pi * r2R * R)) - Vinf2 / 2
        Vitheta2 = ftheta2 / (2 * rho * np.pi * (Vix2 + Vinf))

        if plot == True:
            print(f"Rotor 2 AOA: {np.degrees(alpha)}")
            print(f"Downwash rotor 1: {np.max(Vix)}[m/s]")
            print(f"Downwash rotor 2: {np.max(Vix2)}[m/s")
            print(f"Total Downwash in cruise: {np.max(Vix2+Vix)}[m/s]")
            plt.plot(r2R, funCQi(r2R))
            plt.plot(r2R, funCQ0(r2R))
            # plt.plot(r2R, funcPratio(r2R))
            plt.show()

        return T, pow, Torque

    powmin_1 = 1000000000
    Tmin_1 = 0
    Tmax_1 = 0
    Tmin_2 = 0
    Tmax_2 = 0
    powmax_1 = 0
    powmin_2 = 1000000000
    powmax_2 = 0
    xmin_1 = [0, 0]
    xmax_1 = [0, 0]
    xmin_2 = [0, 0]
    xmax_2 = [0, 0]
    Torquemin_1 = 0
    Torquemin_2 = 0
    Torquemax_1 = 0
    Torquemax_2 = 0

    x0_values = np.linspace(0, 180, 100)
    x1_values = np.linspace(0, np.pi / 2, 100)
    x2_values = np.linspace(0, 180, 100)
    x3_values = np.linspace(0, np.pi / 2, 100)

    Tmax_1 = Tmin_1 = Torquemax_1 = Torquemin_1 = 0
    powmax_1 = powmin_1 = 1e10
    Tmax_2 = Tmin_2 = Torquemax_2 = Torquemin_2 = 0
    powmax_2 = powmin_2 = 1e10

    for x0 in np.linspace(1, 25, 50):
        for x1 in np.linspace(0.01, np.pi / 2, 50):
            for x2 in np.linspace(1, 25, 50):
                for x3 in np.linspace(0.01, np.pi / 2, 50):
                    T_1calc, pow_1calc, Torque_1calc = FirstRotorThrustAndPower(
                        [x0, x1], Vcr, rotorParameters["cruise_AoA"]
                    )
                    T_2calc, pow_2calc, Torque_2calc = SecondRotorThrustAndPower(
                        [x2, x3], T_1calc, Vcr, rotorParameters["cruise_AoA"]
                    )
                    if (
                        np.all(
                            np.array(
                                [T_1calc, pow_1calc, Torque_1calc, T_2calc, pow_2calc, Torque_2calc]
                            )
                            > 0
                        )
                        and (T_1calc + T_2calc) > (2 * Drag)
                        and (pow_1calc + pow_2calc) < (powmin_1)
                        and pow_1calc < 52000
                        and pow_2calc < 115000
                    ):
                        xmin_1 = [x0, x1]
                        xmin_2 = [x2, x3]
                        Tmin_1 = T_1calc
                        Tmin_2 = T_2calc
                        powmin_1 = pow_1calc
                        powmin_2 = pow_2calc
                        Torquemin_1 = Torque_1calc
                        Torquemin_2 = Torque_2calc
                    if (
                        np.all(
                            np.array(
                                [T_1calc, pow_1calc, Torque_1calc, T_2calc, pow_2calc, Torque_2calc]
                            )
                            > 0
                        )
                        and (T_1calc + T_2calc) > (Tmax_1 + Tmax_2)
                        and pow_1calc < 52000
                        and pow_2calc < 115000
                    ):
                        xmax_1 = [x0, x1]
                        xmax_2 = [x2, x3]
                        Tmax_1 = T_1calc
                        Tmax_2 = T_2calc
                        powmax_1 = pow_1calc
                        powmax_2 = pow_2calc
                        Torquemax_1 = Torque_1calc
                        Torquemax_2 = Torque_2calc

    print(f"###################################")
    print(f"First Rotor")
    print(f"Thrust per rotor: {Tmin_1}[N]")
    print(f"Torque per rotor: {Torquemin_1}[Nm]")
    print(f"Power per rotor: {powmin_1/1000} [kW]")
    print(f"Cruise RPM: {xmin_1[0]}")
    print(f"Collective angle: {np.degrees(xmin_1[1])}")
    print(f"MAX CRUISE CONDITIONS")
    print(f"Thrust per rotor: {Tmax_1} [N]")
    print(f"Torque per rotor: {Torquemax_1}[Nm]")
    print(f"Power per rotor: {powmax_1/1000} [kW]")
    print(f"Max RPM: {xmax_1[0]}")
    print(f"Max Collective: {np.degrees(xmax_1[1])}")
    print(f"###################################")
    print(f"Second Rotor")
    print(f"Thrust per rotor: {Tmin_2}[N]")
    print(f"Torque per rotor: {Torquemin_2}[Nm]")
    print(f"Power per rotor: {powmin_2 / 1000} [kW]")
    print(f"Cruise RPM: {xmin_2[0]}")
    print(f"Collective angle: {np.degrees(xmin_2[1])}")
    print(f"MAX CRUISE CONDITIONS")
    print(f"Thrust per rotor: {Tmax_2} [N]")
    print(f"Torque per rotor: {Torquemax_2}[Nm]")
    print(f"Power per rotor: {powmax_2 / 1000} [kW]")
    print(f"Max RPM: {xmax_2[0]}")
    print(f"Max Collective: {np.degrees(xmax_2[1])}")
    print(f"###################################")
    # print(f"Single Rotor")
    # print(f"Thrust per rotor: {TSingleRotor}[N]")
    # print(f"Torque per rotor: {TorqueSingleRotor}[Nm]")
    # print(f"Power per rotor: {powSingleRotor / 1000} [kW]")
    # print(f"Cruise RPM: {xSingleRotor[0] / 10.4 * 2 * np.pi / 60}")
    # print(f"Collective angle: {np.degrees(xSingleRotor[1])}")

    return powmin_1, powmax_1, powmin_2, powmax_2, Tmin_1, Tmax_1, Tmin_2, Tmax_2


def TransitionSimple(V, tilt, cl, S, clWING):
    # SOURCE: PRouty
    R = rotorParameters["Radius"]
    rho = const["airDensity"]
    b = rotorParameters["N_blades"]
    Vmax = const["soundSpeed"] * 0.92
    powmax = rotorParameters["maxEnginePower"]
    Vtang = V * np.cos(tilt)
    omega = (Vmax - Vtang) / R
    mu = V / (omega * R) * np.cos(tilt)
    c = rotorParameters["chord"]
    func = lambda r, phi: (omega * R * (r / R + mu * np.sin(phi))) ** 2 * rho / 2 * cl * c / r
    T = b / (2 * np.pi) * integrate.dblquad(func, 0, 2 * np.pi, 0, R)[0]
    if T * np.sin(tilt) * V > powmax:
        T = powmax / (V * np.sin(tilt))
    L = 1 / 2 * rho * (V + 44 * np.sin(tilt)) ** 2 * S * clWING
    Tvert = T * np.cos(tilt) * 4 + L

    return T * 4, Tvert


def PlotTransition():
    V = np.linspace(0, 120, 100)
    tilt = np.radians(np.linspace(0, 90, 100))
    T = np.empty((100, 100))
    Tvert = np.empty((100, 100))
    i = 0
    for v in V:
        j = 0
        for tlt in tilt:
            Tcalc, Tvertcalc = TransitionSimple(v, tlt, cl=1.25, S=130, clWING=1.63)
            T[j, i] = Tcalc
            Tvert[j, i] = Tvertcalc
            j += 1
        i += 1
    Binary = Tvert >= 2700 * 3.71
    Decel = Tvert <= 2700 * 3.71
    im = plt.pcolormesh(V, np.degrees(tilt), Tvert * Binary, cmap="RdBu")
    clb = plt.colorbar(im)
    clb.ax.set_title("Vertical Thrust [N]")
    plt.xlabel("Velocity")
    plt.ylabel("Tilt Angle")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.show()
    im = plt.pcolormesh(V, np.degrees(tilt), (Tvert - 2700 * 3.71) * Decel / 2700, cmap="RdBu")
    clb = plt.colorbar(im)
    clb.ax.set_title("Acceleration [m/s^2]")
    plt.xlabel("Velocity")
    plt.ylabel("Tilt Angle")
    plt.xlim(left=0, right=120)
    plt.ylim(bottom=0, top=90)
    plt.show()


def AdditionalClimbPower(Drag, V_c):
    MTOW = Drag / 3.71 * 2.205
    R = rotorParameters["Radius"]
    A = np.pi * R**2
    T_hover = Drag
    v_1hover = np.sqrt(T_hover / (2 * const["airDensity"] * A))  # m/s
    v_1climb = -V_c / 2 + np.sqrt((V_c / 2) ** 2 + v_1hover**2)

    Dhp = MTOW / 550 * (V_c / 2 + np.sqrt((V_c / 2) ** 2 + v_1hover**2) - v_1hover)
    collective = (v_1climb + V_c - v_1hover) / (0.75 * (200 / R) * R)

    return Dhp * 745.7, collective


def RotorWash(downwash, H):
    # Source: https://www.tc.faa.gov/its/worldpac/techrpt/rd93-31-1.pdf
    Cd = 1
    rho_rock = 1400
    r_part_max = (downwash**2 * 0.5 * const["airDensity"] * Cd) / (
        rho_rock * const["gravityMars"]
    )

    Beta = (rotorParameters["Radius"] * Cd * const["airDensity"] * 4 * np.pi * r_part_max**2) / (
        4 / 3 * r_part_max**3 * np.pi * rho_rock
    )
    qf = 1 / 2 * const["airDensity"] * downwash**2
    Vparticle = np.sqrt(qf / const["airDensity"]) * (1 - 1 / (np.sqrt(Beta) + 1))
    Eparticle = 0.5 * (4 / 3) * r_part_max**3 * np.pi * rho_rock * Vparticle**2
    print(f"Maximum radius of particle: {r_part_max}[m]")
    print(f"Particle weight: {(4/3)*r_part_max**3*np.pi * rho_rock}[kg]")
    print(f"Maximum particle velocity: {Vparticle} [m/s]")
    print(f"Maximum Energy of Particle:{Eparticle}[J]")

    # particulate cloud

    Kt = 0.025
    R = rotorParameters["Radius"] * 3.28084  # Feet
    DL = 3000 / 4 * 2.20462 / (np.pi * R**2)  # lbs/ft2
    H = H * 3.28084
    rho = const["airDensity"] * 0.00194032  # slugs/ft3
    r2Rj = 2
    r2Rjprev = 1
    while abs((r2Rjprev - r2Rj) / r2Rjprev) > 1e-10:
        t2R = H / R + (r2Rj - 1)
        t2De = 0.707 * t2R
        if t2De < 4:
            qsqn = 1 - 0.025 * t2De**2
        else:
            qsqn = 2.4 / t2De
        Un = (2 * (DL) / rho) ** 0.5
        um = Un * qsqn**0.5
        kg = 1 - 0.9 * np.exp(-2 * H / R)
        U = kg * (Un / 2)
        r2Rjcalc = 2.5081 * (U / um) ** 0.486
        r2Rjprev = r2Rj
        r2Rj = r2Rjcalc

    Um = (0.3586 * (r2Rj) ** 0.885 * um * U**0.14) ** 0.88
    Cu = um / Um * r2Rj**1.143
    Rc = R * (np.sqrt(Kt) / (2.2 * 0.5 * rho * (Um**2) * (Cu**2))) ** (-0.437)

    Rv = 0.785 * Rc
    Zv = 0.329 * Rc
    A = 2 / np.pi * np.log(Zv / (Rc - Rv))
    phi0 = np.pi / 2 * (np.log(Rc - Rv) / (np.log(Zv) - np.log(Rc - Rv)))
    lv = np.exp(A * (-np.pi / 2 + phi0))
    Hcloud = lv + Zv
    print(f"Rc: {Rc/3.28084}[m], Rv: {Rv/3.28084}[m]")
    print(
        f"Rc-Rv: {(Rc-Rv)/3.28084}[m], Radius: {lv/3.28084}[m], Cloud height: {Hcloud/3.28084}[m]"
    )


# DetermineCutout(3000, 200, plot=False)
Trot1, Torquerot1, powrot1, addhprot1 = TakeoffRotor(200, True)
Trot2, Torquerot2, powrot2, addhprot2 = SecondRotor(
    200, rotorParameters["blade_twist"], Trot1, prnt=True
)
RotorWash(44, 10)
# PlotTransition()

powmin_1, powmax_1, powmin_2, powmax_2, Tmin_1, Tmax_1, Tmin_2, Tmax_2 = OptimizeCruise(200, 112)
print(f"########TOTAL AIRCRAFT########")
print(f"Total Lifting Capacity: {2*(Trot1+Trot2)/3.71} [kg]")
print(f"Total Hover Power: {2*(powrot1+powrot2)/1000} [kW]")
print(f"Total Takeoff Power: {2*(powrot1+powrot2+addhprot1+addhprot2)/1000}[kW]")
# print(f'SIngle rotor cruise: {2*powmin_1/1000} [kW] at {2*Tmin_1}')
print(
    f"Total Cruise Power: {2*(powmin_1+powmin_2)/1000} [kW] at a total thrust of {2*(Tmin_1+Tmin_2)}"
)
print(
    f"Total Cruise Power: {2*(powmax_1+powmax_2)/1000} [kW] at a total thrust of {2*(Tmax_1+Tmax_2)}"
)
