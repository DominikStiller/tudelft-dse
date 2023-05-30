import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from scipy.interpolate import InterpolatedUnivariateSpline

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
    "Radius": 10.4,
    'chord': 10.4/20,
    "cutout": 0.1,
    "N_blades": 6,
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
    "blade_twist": 0,
}

def TakeoffRotor(V_tip, prnt, n_elements = 15):
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


    CLalpha = InterpolatedUnivariateSpline(rotorParameters['alpha_data'], rotorParameters['CL_data'])
    CDalpha = InterpolatedUnivariateSpline(rotorParameters['alpha_data'], rotorParameters['CD_data'])
    a0 = CLalpha.derivative(1)

    A = np.pi * R**2

    r2R = np.arange(1, n_elements + 1) / n_elements  # Local radius to total rotor

    c2R = c / R

    M_local = (r2R) * (omega * R / V_m)  # Local Mach number (omega*R = V_tip which is constant)

    a = a0(takeoffAOA) / (np.sqrt(1 - M_local**2))  # Lift curve slope corrected for mach number

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
    rotorParameters["blade_twist"] = theta+rotorParameters["alpha_0"]

    # S1223
    cl = CLalpha(alpha)
    cd = CDalpha(alpha)
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

    Pratiofunc = r2R**3 * (1 - np.sqrt(1 - (2 * Ct / (r2R**2)))) ** 2
    ii = np.isfinite(Pratiofunc)
    funcPratio = InterpolatedUnivariateSpline(r2R[ii], Pratiofunc[ii])

    DCQ_I = 1 / Ct * funcPratio.integral(np.sqrt(2 * Ct), 1) * CQ_induced

    Cq = (CQ_profile + CQ_induced + DCQ_I) / 0.95

    T = rho * A * (omega * R) ** 2 * Ct
    if coaxial:
        T *= 0.88
    pow = rho * A * V_tip**3 * Cq / 550 / 1.341 * 1000
    torque = rho * A * (omega * R) ** 2 * Cq
    # pow = power(T, R)
    sigma = b * c / (np.pi * R)
    ct2sigma = Ct / sigma
    cq2sigma = Cq / sigma


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

    downwash = np.sqrt(T/(2*np.pi*rho*R*R))

    if prnt==True:
        TipVortex(Ct, sigma)
        print(f"Using a rotor radius of: {R}[m]")
        print(f'Cutout: {x0}')
        print(f"Downwash: {downwash} [m/s]")
        print(f"Blade Twist: {np.degrees(rotorParameters['blade_twist'])}")
        print(f"Takeoff Thrust per rotor: {T}[N]")
        print(f"Takeoff Torque per rotor: {torque}[Nm]")
        print(f"Takeoff Power per rotor: {pow / 1000}[kW]")
        print(f"Takeoff RPM: {V_tip/R * 60/(2*np.pi)}")
        print(f"Rotor Mass: {Rotor_mass}")
        print(f"Lifting capacity: {T * 4 / gm}[kg]")

    return T, torque


def ClimbingVertical(MTOM,Vclimb):
    T = MTOM * const["gravityMars"]
    rho = const["airDensity"]
    V_sound = const["soundSpeed"]
    R = rotorParameters["Radius"]
    omega = V_sound / R * 0.92

    v_i_hover = np.sqrt(T / (2 * rho * R**2 * np.pi))

    GW = 2700 * 2.20462
    v_i_climb = np.sqrt((2700 * 3.71) / (2 * rho * R**2 * np.pi))
    Dcollective = (v_i_climb + Vclimb - v_i_hover) / (0.75 * omega * R)

    v_i_hover *= 3.28084
    v_i_climb *= 3.28084
    Vclimb *= 3.28084

    Dhp = GW / 550 * (Vclimb / 2 + np.sqrt((Vclimb / 2) ** 2 + v_i_hover**2) - v_i_hover)

    print(f"Required additional power: {Dhp}[hp]")
    print(f"Required collective angle: {np.degrees(Dcollective)} [deg]")

def TipVortex(Ct, sigma):
    v = const["visc_cr"]
    vt = 9.82e-6
    R = rotorParameters["Radius"]
    c = rotorParameters["chord"]
    v_tip = const["soundSpeed"]*0.92
    omega = v_tip/R
    Re_vortex = 2*omega*R*c/v * Ct/sigma

    r0 = 0.029 #deg
    alphaL = 1.25643 #lamb constant
    delta = 1+ vt/v

    t = np.linspace(0, 5, 100)

    rc = np.sqrt(r0**2+4*alphaL*v*delta*t)

    #plt.plot(t, rc)
    #plt.show()

    print(f"Vortex tip Reynolds Number: {Re_vortex}")

def RotorThrust(N_blades, R, V_tip, coaxial, cutout):
    # Rename constants
    gm = const["gravityMars"]
    rho = const["airDensity"]
    V_m = const["soundSpeed"]

    alpha = rotorParameters["alpha_data"]
    CL = rotorParameters["CL_data"]
    CD = rotorParameters["CD_data"]

    CL2alpha = InterpolatedUnivariateSpline(alpha, CL)
    CD2alpha = InterpolatedUnivariateSpline(alpha, CD)

    # if N_rotors <= 0:
    #     return "N_rotors has to be greater than zero.",0,0,0
    # elif N_blades <= 0:
    #     return "N_blades has to be greater than zero.",0,0,0

    b = N_blades
    v_tip = V_tip

    c = R / 20  # Update the chord

    x0 = cutout  # Distance from the blade's base to the centre of rotation
    omega = v_tip / R * 0.92  # Angular velocity

    n_elements = 15  # Split up blade into set of elements
    a0 = 6
    alpha0 = -np.radians(5)
    A = np.pi * R**2

    r2R = np.arange(1, n_elements + 1) / n_elements  # Local radius to total rotor

    c2R = c / R

    M_local = (r2R) * (omega * R / V_m)  # Local Mach number (omega*R = V_tip which is constant)

    a = a0 / (np.sqrt(1 - M_local**2))  # Lift curve slope corrected for mach number

    def v12Omegar(theta):
        return (
            (a * b * c2R)
            / (16 * np.pi * r2R)
            * (-1 + np.sqrt(1 + (32 * np.pi * theta * r2R) / (a * b * c2R)))
        )

    alpha = np.radians(6) * np.ones(n_elements)

    def func(x):
        return alpha - (x - np.arctan(v12Omegar(x)))

    theta = scipy.optimize.fsolve(
        func,
        0.1
        * np.ones(
            n_elements,
        ),
    )
    Dcollective = np.radians(np.linspace(-15, 5, 100))
    T_list = np.empty(np.shape(Dcollective))
    Torque_list = np.empty(np.shape(Dcollective))
    i = 0
    for Dtheta in Dcollective:
        theta_local = theta + Dtheta
        aoa = theta_local - np.arctan(v12Omegar(theta_local))
        cl_loc = CL2alpha(np.degrees(aoa))
        cd_loc = CD2alpha(np.degrees(aoa))

        DctDr2R = b * r2R**2 * c2R * cl_loc / (2 * np.pi)
        # Create spline of data points in order
        funCT = InterpolatedUnivariateSpline(r2R, DctDr2R)
        Ct_lossless = funCT.integral(x0, 1)
        # Loss of lift from the tip of the rotor
        if Ct_lossless < 0.006:
            B = 1 - 0.06 / b
        else:
            B = 1 - np.sqrt(2.27 * Ct_lossless - 0.01) / b

        Ct = Ct_lossless - funCT.integral(B, 1)

        Dcq0Dr2R = b * r2R ** 3 * c2R * cd_loc / (2 * np.pi)
        funCQ0 = InterpolatedUnivariateSpline(r2R, Dcq0Dr2R)
        CQ_profile = funCQ0.integral(x0, 1)
        DcqiDr2R = b * r2R ** 3 * c2R * cl_loc * v12Omegar(theta) / (2 * np.pi)
        funCQi = InterpolatedUnivariateSpline(r2R, DcqiDr2R)
        CQ_induced = funCQi.integral(B, 1)

        Pratiofunc = r2R ** 3 * (1 - np.sqrt(1 - (2 * Ct / (r2R ** 2)))) ** 2
        ii = np.isfinite(Pratiofunc)
        if any(ii)==True:
            funcPratio = InterpolatedUnivariateSpline(r2R[ii], Pratiofunc[ii])

            DCQ_I = 1 / Ct * funcPratio.integral(np.sqrt(2 * Ct), 1) * CQ_induced
        else:
            DCQ_I=0

        Cq = (CQ_profile + CQ_induced + DCQ_I) / 0.95

        T_list[i] = rho * A * (omega * R) ** 2 * Ct
        if coaxial:
            T_list[i] *= 0.88
        Torque_list[i] = rho * A * (omega * R) ** 2 * Cq
        i+=1
    return V_tip * np.ones(np.shape(Dcollective)), Dcollective, T_list, Torque_list

def RPMCollective(x0):
    V = []
    theta = []
    thrust = []
    torque = []
    for V_tip in np.linspace(10, 220, 210):
        v, th, thrst, trq = RotorThrust(6, 10.4, V_tip, coaxial=True, cutout=x0)
        V = np.append(V, v)
        theta = np.append(theta, th)
        thrust = np.append(thrust, thrst)
        torque = np.append(torque, trq)
    plt.scatter((V/10.4)*9.5492968, np.degrees(theta), c=thrust, cmap=plt.get_cmap("plasma"))
    plt.xlabel('Rotation speed [RPM]')
    plt.ylabel('Collective Angle [deg]')
    plt.colorbar(label="Thrust [N]")
    plt.show()
    plt.scatter((V/10.4)*9.5492968, np.degrees(theta), c=torque, cmap=plt.get_cmap("plasma"))
    plt.xlabel('Rotation speed [RPM]')
    plt.ylabel('Collective Angle [deg]')
    plt.colorbar(label="Torque [Nm]")
    plt.show()

def DetermineCutout(MTOM, tipspeed, plot=False):
    T = 4000*np.ones(100)
    Torque = np.zeros(100)
    i = 0

    while not min(T)<=(MTOM*3.71/4):
        rotorParameters["cutout"]+=0.01
        T[i], Torque[i] = TakeoffRotor(tipspeed, prnt=False, n_elements=15)
        i += 1
    T=T[:i]
    Torque = Torque[:i]
    if plot==True:
        plt.scatter(np.linspace(0, rotorParameters["cutout"], len(T)), T)
        plt.show()
        plt.scatter(np.linspace(0, rotorParameters["cutout"], len(T)), Torque)
        plt.show()
    print(f"Max possible cutout to lift a mass of {MTOM}[kg] is {rotorParameters['cutout']*100}%")

def RotorPerfCruise(n_elements, DragperRotor, prnt):
    # Rename constants
    gm = const["gravityMars"]
    rho = const["airDensity"]
    V_m = const["soundSpeed"]

    b = rotorParameters["N_blades"]
    R= rotorParameters["Radius"]
    c= rotorParameters["chord"]
    x0 = rotorParameters["cutout"]
    aoa = rotorParameters["cruise_AoA"]
    aoa = np.radians(-4.5)
    theta0 = rotorParameters["blade_twist"]
    coaxial = rotorParameters["coaxial"]


    alpha_data = rotorParameters["alpha_data"]
    CL = rotorParameters["CL_data"]

    CD = rotorParameters["CD_data"]

    CLalpha = InterpolatedUnivariateSpline(alpha_data, CL)
    CDalpha = InterpolatedUnivariateSpline(alpha_data, CD)
    a0 = CLalpha.derivative(1)

    A = np.pi * R**2

    r2R = np.arange(1, n_elements + 1) / n_elements  # Local radius to total rotor

    c2R = c / R
    v_tip = 10 #m/s
    T=0
    while T<DragperRotor:
        v_tip+=0.1
        omega = v_tip / R
        M_local = (r2R) * (omega * R / V_m)  # Local Mach number (omega*R = V_tip which is constant)

        a = a0(aoa) / (np.sqrt(1 - M_local**2))  # Lift curve slope corrected for mach number

        def v12Omegar(theta):
            V_c=112
            v_1hov = (
                (a * b * c2R)
                / (16 * np.pi * r2R)
                * (-1 + np.sqrt(1 + (32 * np.pi * theta * r2R) / (a * b * c2R)))
            )

            return ((-omega / 2 * a * c * b + 4 * np.pi * V_c) + np.sqrt(
        (omega / 2 * a * c * b + 4 * np.pi * V_c) ** 2 + 8 * np.pi * b * omega ** 2 * a * c * r2R * (
                    theta - V_c / (omega * r2R))))/(8*np.pi*omega*r2R)
        def func(x):
            theta=theta0+x-rotorParameters['alpha_0']
            return aoa-(theta-np.arctan(v12Omegar(theta)))
        collective = scipy.optimize.fsolve(func, np.ones(15)*0.1)

        theta = theta0+np.average(collective)-rotorParameters['alpha_0']
        alpha = theta-np.arctan(v12Omegar(theta))


        # S1223
        cl = CLalpha(alpha)
        cd = CDalpha(alpha)
        DctDr2R = b * r2R**2 * c2R * cl / (2 * np.pi)
        # Create spline of data points in order
        funCT = InterpolatedUnivariateSpline(r2R[1:], DctDr2R[1:])
        Ct_lossless = funCT.integral(x0, 1)
        # Loss of lift from the tip of the rotor
        if Ct_lossless < 0.006:
            B = 1 - 0.06 / b
        else:
            B = 1 - np.sqrt(2.27 * Ct_lossless - 0.01) / b

        Ct = Ct_lossless - funCT.integral(B, 1)

        Dcq0Dr2R = b * r2R**3 * c2R * cd / (2 * np.pi)
        funCQ0 = InterpolatedUnivariateSpline(r2R[1:], Dcq0Dr2R[1:])
        CQ_profile = funCQ0.integral(x0, 1)
        DcqiDr2R = b * r2R**3 * c2R * cl * v12Omegar(theta) / (2 * np.pi)

        funCQi = InterpolatedUnivariateSpline(r2R[1:], DcqiDr2R[1:])
        CQ_induced = funCQi.integral(B, 1)

        Pratiofunc = r2R**3 * (1 - np.sqrt(1 - (2 * Ct / (r2R**2)))) ** 2
        ii = np.isfinite(Pratiofunc)
        funcPratio = InterpolatedUnivariateSpline(r2R[ii], Pratiofunc[ii])

        DCQ_I = 1 / Ct * funcPratio.integral(np.sqrt(2 * Ct), 1) * CQ_induced

        Cq = (CQ_profile + CQ_induced + DCQ_I) / 0.95

        T = rho * A * (omega * R) ** 2 * Ct
        if coaxial:
            T *= 0.88
        pow = rho * A * v_tip**3 * Cq / 550 / 1.341 * 1000

        torque = rho * A * (omega * R) ** 2 * Cq
        sigma = b * c / (np.pi * R)
        ct2sigma = Ct / sigma
        cq2sigma = Cq / sigma
    pow += AdditionalClimbPower(T, const["cruiseSpeed"])
    RPM = (v_tip/10.4)*60/(2*np.pi)

    plt.plot(r2R,theta)
    plt.show()

    downwash = np.sqrt(T/(2*np.pi*rho*R**2))

    V_c = 112

    v_1 = ((-omega / 2 * a * c * b + 4 * np.pi * V_c) + np.sqrt(
        (omega / 2 * a * c * b + 4 * np.pi * V_c) ** 2 + 8 * np.pi * b * omega ** 2 * a * c * r2R * (
                    theta - V_c / (omega * r2R))))/(8*np.pi)
    print(v_1)
    if prnt==True:
        print(f"Cruise Collective angle: {np.degrees(np.mean(collective))}")
        print(f"Downwash: {downwash}[m/s]")
        print(f"Cruise Thrust per rotor: {T}[N]")
        print(f"Cruise Torque per rotor: {torque}[Nm]")
        print(f"Cruise Power per rotor: {pow / 1000} or {T*const['cruiseSpeed']/1000}[kW]")
        print(f"Cruise RPM: {RPM}")

    return T, torque

def AdditionalClimbPower(Drag, V_c):
    MTOW = Drag/3.71* 2.205
    R = rotorParameters["Radius"]
    A = np.pi * R**2
    T_max = 3000*3.71 / 4
    T_hover = 800 / 4
    v_1hover = np.sqrt(T_hover / (2 * const["airDensity"] * A))  # m/s

    Dhp = MTOW / 550 * (V_c / 2 + np.sqrt((V_c / 2) ** 2 + v_1hover**2) - v_1hover)
    return Dhp*745.7


DetermineCutout(6600, 210)
TakeoffRotor(210, True)
RotorPerfCruise(15, DragperRotor=200, prnt=True)