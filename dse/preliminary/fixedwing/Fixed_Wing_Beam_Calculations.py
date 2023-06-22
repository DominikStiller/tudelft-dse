import csv

import matplotlib.pyplot as plt
import numpy as np

from dse.plotting import format_plot, save_plot

filenameArray = [
    "Flying_Wing_Data_AoA=9_3.csv",
    "Bi_Wing_1_a=14.00_v=111.00ms.csv",
    "Bi_Wing_2_a=14.00_v=111.00ms.csv",
    "S1223.csv",
]
n = 10000

### Cruise conditions
V_cruise = 400 / 3.6  # [m/s]
rho = 0.01  # [kg/m3] air density
q = 0.5 * rho * (V_cruise**2)  # [Pa] Dynamic pressure
n = 2.5  # [-] Maximum load factor


def GetAirfoil(filename):
    ### Airfoil Data
    with open(filename) as csvfile:
        data = csv.reader(csvfile, delimiter="\t")
        Data = []
        for row in data:
            Data.append(row)
    Airfoilcoordinates = [[float(Data[i][0]), float(Data[i][1])] for i in range(len(Data))]
    return Airfoilcoordinates


def GetWingbox(Indx):  # 0 and 1 for airfoil wingbox, 2 for square wingbox
    WingboxArray = np.zeros((6, 4))
    if Indx == 2:
        # for cantilever beam
        WingboxArray[0] = np.array([0.0, 1, 1, 0.0])  # [-] As percentage of unit chord
        WingboxArray[1] = np.array([0.5, 0.5, -0.5, -0.5])  # [-] As percentage of unit chord
        WingboxArray[2] = np.array([0.001, 0.001, 0.001, 0.001])
    else:
        WingboxArray[0] = np.array([0.05, 0.5, 0.5, 0.05])  # [-] As percentage of unit chord
        WingboxArray[1] = np.array(
            [0.07, 0.12486, 0.04953, -0.01123]
        )  # [-] As percentage of unit chord
        WingboxArray[2] = np.array([0.001, 0.01, 0.001, 0.01])  # [m] [t12, t23, t34, t41
    b_lst = []  # [-] The length of each sheet for the unit chord: [b12, b23, b34, b41]
    for i in range(len(WingboxArray[0]) - 1):
        b_lst.append(
            (
                (WingboxArray[0, i] - WingboxArray[0, i + 1]) ** 2
                + (WingboxArray[1, i] - WingboxArray[1, i + 1]) ** 2
            )
            ** 0.5
        )
    b_lst.append(
        (
            (WingboxArray[0, 0] - WingboxArray[0, -1]) ** 2
            + (WingboxArray[1, 0] - WingboxArray[1, -1]) ** 2
        )
        ** 0.5
    )
    WingboxArray[3] = np.array(b_lst)
    xcg_lst = []
    ycg_lst = []
    for i in range(len(WingboxArray[0]) - 1):
        xcg_lst.append(WingboxArray[0, i] + (WingboxArray[0, i + 1] - WingboxArray[0, i]) / 2)
        ycg_lst.append(WingboxArray[1, i] + (WingboxArray[1, i + 1] - WingboxArray[1, i]) / 2)
    xcg_lst.append(WingboxArray[0, -1] + (WingboxArray[0, 0] - WingboxArray[0, -1]) / 2)
    ycg_lst.append(WingboxArray[1, -1] + (WingboxArray[1, 0] - WingboxArray[1, -1]) / 2)
    WingboxArray[4] = np.array(xcg_lst)
    WingboxArray[5] = np.array(ycg_lst)
    return WingboxArray


def NeutralAxisUnit(WingboxArray):
    xbar = sum(
        WingboxArray[3, i] * WingboxArray[2, i] * WingboxArray[4, i]
        for i in range(len(WingboxArray[3]))
    ) / sum(WingboxArray[3, i] * WingboxArray[2, i] for i in range(len(WingboxArray[3])))
    ybar = sum(
        WingboxArray[3, i] * WingboxArray[2, i] * WingboxArray[5, i]
        for i in range(len(WingboxArray[3]))
    ) / sum(WingboxArray[3, i] * WingboxArray[2, i] for i in range(len(WingboxArray[3])))
    NA = np.array([xbar, ybar])
    return NA


def IxxCalcUnit(WingboxArray, NA):  # Multiply by c^3 to gain the moment of inertia at said location
    Ixx12 = (
        WingboxArray[2, 0]
        * WingboxArray[3, 0] ** 3
        * ((WingboxArray[1, 1] - WingboxArray[1, 0]) / WingboxArray[3, 0]) ** 2
    ) / 12
    Ixx23 = (WingboxArray[2, 1] * WingboxArray[3, 1] ** 3) / 12
    Ixx34 = (
        WingboxArray[2, 2]
        * WingboxArray[3, 2] ** 3
        * ((WingboxArray[1, 3] - WingboxArray[1, 2]) / WingboxArray[3, 2]) ** 2
    ) / 12
    Ixx41 = (WingboxArray[2, 3] * WingboxArray[3, 3] ** 3) / 12
    IxxPart = [Ixx12, Ixx23, Ixx34, Ixx41]
    Ixx = sum(
        IxxPart[i] + WingboxArray[2, i] * WingboxArray[3, i] * (WingboxArray[5, i] - NA[1]) ** 2
        for i in range(len(IxxPart))
    )
    print(Ixx)
    return Ixx


def IyyCalcUnit(WingboxArray, NA):  # Multiply by c^4 to gain the moment of inertia at said location
    Iyy12 = (
        WingboxArray[2, 0]
        * WingboxArray[3, 0] ** 3
        * ((WingboxArray[0, 1] - WingboxArray[0, 0]) / WingboxArray[3, 0]) ** 2
    ) / 12
    Iyy23 = 0
    Iyy34 = (
        WingboxArray[2, 2]
        * WingboxArray[3, 2] ** 3
        * ((WingboxArray[0, 3] - WingboxArray[0, 2]) / WingboxArray[3, 2]) ** 2
    ) / 12
    Iyy41 = 0
    IyyPart = [Iyy12, Iyy23, Iyy34, Iyy41]
    Iyy = sum(
        IyyPart[i] + WingboxArray[2, i] * WingboxArray[3, i] * (WingboxArray[4, i] - NA[0]) ** 2
        for i in range(len(IyyPart))
    )
    return Iyy


def IxyCalcUnit(WingboxArray, NA):  # Multiply by c^4 to gain the moment of inertia at said location
    Ixy12 = (
        WingboxArray[2, 0]
        * WingboxArray[3, 0] ** 3
        * (
            ((WingboxArray[1, 1] - WingboxArray[1, 0]) * (WingboxArray[0, 1] - WingboxArray[0, 0]))
            / (WingboxArray[3, 0] ** 2)
        )
    ) / 12
    Ixy23 = 0
    Ixy34 = (
        WingboxArray[2, 2]
        * WingboxArray[3, 2] ** 3
        * (
            ((WingboxArray[1, 3] - WingboxArray[1, 2]) * (WingboxArray[0, 3] - WingboxArray[0, 2]))
            / (WingboxArray[3, 2] ** 2)
        )
    ) / 12
    Ixy41 = 0
    IxyPart = [Ixy12, Ixy23, Ixy34, Ixy41]
    Ixy = sum(
        IxyPart[i]
        + WingboxArray[2, i]
        * WingboxArray[3, i]
        * (WingboxArray[5, i] - NA[1])
        * (WingboxArray[4, i] - NA[0])
        for i in range(len(IxyPart))
    )
    return Ixy


def GetDataXFLR(filename):
    with open(filename) as csvfile:
        data = csv.reader(csvfile, delimiter=",")
        Data = []
        for row in data:
            Data.append(row)

    y_lst = [float(i[0]) for i in Data[1:]]  # [m] Y-step of the applied loads b
    c_lst = [float(i[1]) for i in Data[1:]]  # [m] Y-step of the applied loads b
    AoA_lst = [float(i[2]) for i in Data[1:]]  # [degrees] induced angle of angle of attack
    CL_lst = [float(i[3]) for i in Data[1:]]  # [-] Lift coefficient of the wing
    PCd_lst = [float(i[4]) for i in Data[1:]]  # [-] Parasite drag coefficient of the wing
    ICd_lst = [float(i[5]) for i in Data[1:]]  # [-] Induced drag coefficient of the wing
    CmGeom_lst = [float(i[6]) for i in Data[1:]]  # [-] Moment coefficient
    CmAirfchord4_lst = [float(i[7]) for i in Data[1:]]  # [-] Moment coefficient at quarter chord
    Data = np.zeros((8, len(y_lst)))
    Data[0] = np.array(y_lst)
    Data[1] = np.array(c_lst)
    Data[2] = np.array(AoA_lst)
    Data[3] = np.array(CL_lst)
    Data[4] = np.array(PCd_lst)
    Data[5] = np.array(ICd_lst)
    Data[6] = np.array(CmGeom_lst)
    Data[7] = np.array(CmAirfchord4_lst)
    return Data


def ConstructDataCB(n):
    y_lst = np.linspace(-19.025, -0.4982, n)
    c = 1
    c_lst = np.ones(len(y_lst)) * c
    PCd_lst = np.zeros(len(y_lst))
    ICd_lst = np.zeros(len(y_lst))
    Data = np.zeros((8, len(y_lst)))
    Data[0] = np.array(y_lst)

    Data[1] = np.array(c_lst)
    Data[2] = np.zeros(n)
    Data[3] = np.zeros(n)
    Data[4] = np.array(PCd_lst)
    Data[5] = np.array(ICd_lst)
    Data[6] = np.zeros(n)
    Data[7] = np.zeros(n)
    return Data


def AppliedLoadsCB(Data):
    Loads = np.ones((6, len(Data[0])))
    Mx_lst = [0]
    My_lst = np.zeros(len(Data[0]))
    Mz_lst = np.zeros(len(Data[0]))
    L_lst = [0]
    D_lst = np.zeros(len(Data[0]))
    w_lst = []
    Y_lst = [0]
    for i in range(len(Data[0]) - 1):
        dy = Data[0, i + 1] - Data[0, i]
        L = 300
        w = 300 * dy
        w_lst.append(w)
        L_lst.append(L)
        Mx_lst.append(Mx_lst[i] + w * dy / 2 + sum(w_lst) * dy)
        Y_lst.append(1)  # place holder for axial forces for rotor blades.

    Loads[0] = np.array(Mx_lst)
    Loads[1] = np.array(My_lst)
    Loads[2] = np.array(Mz_lst)
    Loads[3] = np.array(L_lst)
    Loads[4] = np.array(D_lst)
    Loads[5] = np.array(Y_lst)
    return Loads


def AppliedLoads(Data):
    Loads = np.ones((6, len(Data[0])))
    Mx_lst = [0]
    My_lst = [0]
    Mz_lst = [0]
    L_lst = [0]
    D_lst = [0]
    Y_lst = [0]
    for i in range(len(Data[0]) - 1):
        dy = Data[0, i + 1] - Data[0, i]
        c_avg = (Data[1, i + 1] + Data[1, i]) / 2
        L = q * Data[3, i] * c_avg * dy
        D = q * (Data[5, i] + Data[4, i]) * c_avg * dy
        Y_lst.append(1)  # place holder for axial forces for rotor blades.
        L_lst.append(L)
        D_lst.append(D)
        Mx_lst.append(Mx_lst[i] + L * dy / 2 + sum(L_lst) * dy)
        My_lst.append(My_lst[i] + D * dy / 2 + sum(D_lst) * dy)
        Mz_lst.append(1)  # Apply torsion for thinwalled structure
    Loads[0] = np.array(Mx_lst)
    Loads[1] = np.array(My_lst)
    Loads[2] = np.array(Mz_lst)
    Loads[3] = np.array(L_lst)
    Loads[4] = np.array(D_lst)
    Loads[5] = np.array(Y_lst)
    return Loads


def InternalStress(Loads, WingboxArray, Data, NA):
    Ixx = IxxCalcUnit(WingboxArray, NA)
    Iyy = IyyCalcUnit(WingboxArray, NA)
    Ixy = IxyCalcUnit(WingboxArray, NA)
    Stress = np.ones((6, 4, len(Data[0])))
    sigmax = np.ones((len(WingboxArray[0]), len(Data[0])))  # Add formulas for these
    sigmay = np.ones((len(WingboxArray[0]), len(Data[0])))
    sigmaz = np.ones((len(WingboxArray[0]), len(Data[0])))  # Add formulas for these
    tauxy = np.ones((len(WingboxArray[0]), len(Data[0])))  # Add formulas for these
    tauyz = np.ones((len(WingboxArray[0]), len(Data[0])))  # Add formulas for these
    tauxz = np.ones((len(WingboxArray[0]), len(Data[0])))  # Add formulas for these
    sigmay_1_lst = [
        10**-6
        * (
            (Loads[0, i] * Iyy - Loads[1, i] * Ixy) * (-WingboxArray[1, 0] + NA[1])
            + (Loads[1, i] * Ixx - Loads[0, i] * Ixy) * -(-WingboxArray[0, 0] + NA[0])
        )
        / ((Ixx * Iyy - Ixy**2) * Data[1, i] ** 2)
        for i in range(len(Data[0]))
    ]
    sigmay_2_lst = [
        10**-6
        * (
            (Loads[0, i] * Iyy - Loads[1, i] * Ixy) * (-WingboxArray[1, 1] + NA[1])
            + (Loads[1, i] * Ixx - Loads[0, i] * Ixy) * (WingboxArray[0, 1] - NA[0])
        )
        / ((Ixx * Iyy - Ixy**2) * Data[1, i] ** 2)
        for i in range(len(Data[0]))
    ]
    sigmay_3_lst = [
        10**-6
        * (
            (Loads[0, i] * Iyy - Loads[1, i] * Ixy) * (-WingboxArray[1, 2] + NA[1])
            + (Loads[1, i] * Ixx - Loads[0, i] * Ixy) * -(-WingboxArray[0, 2] + NA[0])
        )
        / ((Ixx * Iyy - Ixy**2) * Data[1, i] ** 2)
        for i in range(len(Data[0]))
    ]
    sigmay_4_lst = [
        10**-6
        * (
            (Loads[0, i] * Iyy - Loads[1, i] * Ixy) * (-WingboxArray[1, 3] + NA[1])
            + (Loads[1, i] * Ixx - Loads[0, i] * Ixy) * (WingboxArray[0, 3] - NA[0])
        )
        / ((Ixx * Iyy - Ixy**2) * Data[1, i] ** 2)
        for i in range(len(Data[0]))
    ]
    sigmay[0] = np.array(sigmay_1_lst)
    sigmay[1] = np.array(sigmay_2_lst)
    sigmay[2] = np.array(sigmay_3_lst)
    sigmay[3] = np.array(sigmay_4_lst)
    Stress[0] = sigmax
    Stress[1] = sigmay
    Stress[2] = sigmaz
    Stress[3] = tauxy
    Stress[4] = tauyz
    Stress[5] = tauxz
    return Stress


def analyticalSol(Data):
    b2 = Data[0, 0]
    b = 1
    h = 1
    t = 0.001
    Ixx = ((b + t) * (h + t) ** 3 - (b - t) * (h - t) ** 3) / 12
    w = 300
    Mx = w / 2 * (b2 - Data[0]) ** 2
    z_booms = np.ones((4, np.size(Data[0]))) * 0.5
    z_booms[2:] *= -1.0
    sigmaA = Mx * z_booms / Ixx * 10**-6
    print("k")
    print(Ixx)
    return Mx, sigmaA


def MakePlots0(indx, WingboxArray, Data, Loads, NA):
    fig, ax1 = plt.subplots(figsize=(8, 4))

    color = "tab:red"
    ax1.set_xlabel("Span [m]")
    ax1.set_ylabel("Lift [N]", color=color)
    ax1.plot(-Data[0], Loads[3], color=color)
    ax1.tick_params(axis="y", labelcolor=color)

    ax2 = ax1.twinx()

    color = "tab:blue"
    ax2.set_ylabel("Drag [N]", color=color)  # we already handled the x-label with ax1
    ax2.plot(-Data[0], Loads[4], color=color, linestyle="--")
    ax2.tick_params(axis="y", labelcolor=color)

    format_plot()
    save_plot(".", "AppliedLoadsFlyingWing")
    plt.show()

    # Making the internal loads plots

    fig, ax1 = plt.subplots(figsize=(8, 4))

    color = "tab:red"
    ax1.set_xlabel("Span [m]")
    ax1.set_ylabel("Mx [Nm]", color=color)
    ax1.plot(-Data[0], Loads[0], color=color)
    ax1.tick_params(axis="y", labelcolor=color)

    ax2 = ax1.twinx()

    color = "tab:blue"
    ax2.set_ylabel("Mz [Nm]", color=color)  # we already handled the x-label with ax1
    ax2.plot(-Data[0], Loads[1], color=color, linestyle="--")
    ax2.tick_params(axis="y", labelcolor=color)

    format_plot()
    save_plot(".", "InternalLoadsFlyingwing")
    plt.show()

    # # Internal stress diagram
    # sigmaz_1_lst, sigmaz_2_lst, sigmaz_3_lst, sigmaz_4_lst = InternalStress()
    plt.figure(figsize=(8, 4))
    plt.plot(-Data[0], Stress[1, 0], label="Stress distribution in point 1")
    plt.plot(-Data[0], Stress[1, 1], label="Stress distribution in point 2")
    plt.plot(-Data[0], Stress[1, 2], label="Stress distribution in point 3")
    plt.plot(-Data[0], Stress[1, 3], label="Stress distribution in point 4")
    plt.legend()
    plt.ylabel("Stress [Mpa]")
    plt.xlabel("Span [m]")
    format_plot()
    save_plot(".", "StressesFlyingWing")
    plt.show()

    # Plot the wingbox structure
    fig, ax = plt.subplots(figsize=(8, 4))
    # xbar, NA[1] = NeutralAxisUnit()

    ax.axhline(NA[1], linestyle="--", color="darkgrey", zorder=1)
    ax.axvline(NA[0], linestyle="--", color="darkgrey", zorder=2)
    if indx != 2:
        ax.plot(*list(zip(*Airfoilcoordinates)), color="darkslategrey", zorder=3)
    ax.plot(WingboxArray[0], WingboxArray[1], color="black", zorder=4)
    ax.plot(
        [WingboxArray[0, 0], WingboxArray[0, -1]],
        [WingboxArray[1, 0], WingboxArray[1, -1]],
        color="black",
        zorder=4,
    )
    n = [1, 2, 3, 4]

    ax.scatter(WingboxArray[0], WingboxArray[1], color="orange", zorder=5)
    for i, txt in enumerate(n):
        ax.annotate(txt, (WingboxArray[0, i], WingboxArray[1, i]), zorder=6)
    plt.gca().axis("equal")
    plt.ylabel("z/MAC [-]")
    plt.xlabel("x/MAC [-]")
    format_plot()
    save_plot(".", "WingboxStructure")
    plt.show()


def MakePlots1(indx, WingboxArray, Data1, Data2, Loads1, loads2, NA):
    fig, ax1 = plt.subplots(figsize=(8, 4))

    color = "tab:red"
    ax1.set_xlabel("Span [m]")
    ax1.set_ylabel("Lift [N]", color=color)
    ax1.plot(-Data1[0], Loads1[3], color=color, label="Upper wing")
    ax1.plot(-Data2[0], Loads2[3], color=color, linestyle="--", label="Lower wing")
    ax1.tick_params(axis="y", labelcolor=color)

    ax2 = ax1.twinx()

    color = "tab:blue"
    ax2.set_ylabel("Drag [N]", color=color)  # we already handled the x-label with ax1
    ax2.plot(-Data1[0], Loads1[4], color=color, label="Upper wing")
    ax2.plot(-Data2[0], Loads2[4], color=color, linestyle="--", label="Lower wing")
    ax2.tick_params(axis="y", labelcolor=color)

    format_plot()
    save_plot(".", "AppliedLoadsBiplane")
    plt.show()

    # Making the internal loads plots

    fig, ax1 = plt.subplots(figsize=(8, 4))

    color = "tab:red"
    ax1.set_xlabel("Span [m]")
    ax1.set_ylabel("Mx [Nm]", color=color)
    ax1.plot(-Data1[0], Loads1[0], color=color, label="Upper wing")
    ax1.plot(-Data2[0], Loads2[0], color=color, linestyle="--", label="Lower wing")
    ax1.tick_params(axis="y", labelcolor=color)

    ax2 = ax1.twinx()

    color = "tab:blue"
    ax2.set_ylabel("Mz [Nm]", color=color)  # we already handled the x-label with ax1
    ax2.plot(-Data1[0], Loads1[1], color=color, label="Upper wing")
    ax2.plot(-Data2[0], Loads2[1], color=color, linestyle="--", label="Lower wing")
    ax2.tick_params(axis="y", labelcolor=color)

    format_plot()
    save_plot(".", "InternalLoadsBiplane")
    plt.show()

    # # Internal stress diagram
    # sigmaz_1_lst, sigmaz_2_lst, sigmaz_3_lst, sigmaz_4_lst = InternalStress()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
    ax1.set_title("Upper wing")
    ax1.plot(-Data1[0], Stress1[1, 0], label="Stress distribution in point 1")
    ax1.plot(-Data1[0], Stress1[1, 1], label="Stress distribution in point 2")
    ax1.plot(-Data1[0], Stress1[1, 2], label="Stress distribution in point 3")
    ax1.plot(-Data1[0], Stress1[1, 3], label="Stress distribution in point 4")
    ax1.legend()

    ax1.set_ylabel("Stress [Mpa]")
    ax2.set_title("Lower wing")
    ax2.plot(-Data2[0], Stress2[1, 0], label="Stress distribution in point 1")
    ax2.plot(-Data2[0], Stress2[1, 1], label="Stress distribution in point 2")
    ax2.plot(-Data2[0], Stress2[1, 2], label="Stress distribution in point 3")
    ax2.plot(-Data2[0], Stress2[1, 3], label="Stress distribution in point 4")
    ax2.legend()

    ax2.set_ylabel("Stress [Mpa]")
    ax2.set_xlabel("Span [m]")
    format_plot()
    save_plot(".", "StressesBiwing")
    plt.show()

    # Plot the wingbox structure
    fig, ax = plt.subplots()
    # xbar, NA[1] = NeutralAxisUnit()

    ax.axhline(NA[1], linestyle="--", color="darkgrey", zorder=1)
    ax.axvline(NA[0], linestyle="--", color="darkgrey", zorder=2)
    if indx != 2:
        ax.plot(*list(zip(*Airfoilcoordinates)), color="lightgrey", zorder=3)
    ax.plot(WingboxArray[0], WingboxArray[1], color="black", zorder=4)
    ax.plot(
        [WingboxArray[0, 0], WingboxArray[0, -1]],
        [WingboxArray[1, 0], WingboxArray[1, -1]],
        color="black",
        zorder=4,
    )
    n = [1, 2, 3, 4]

    ax.scatter(WingboxArray[0], WingboxArray[1], color="orange", zorder=5)
    for i, txt in enumerate(n):
        ax.annotate(txt, (WingboxArray[0, i], WingboxArray[1, i]), zorder=6)
    plt.gca().axis("equal")
    plt.ylabel("z/MAC [-]")
    plt.xlabel("x/MAC [-]")
    plt.show()


def MakePlots2(indx, WingboxArray, Data, Loads, NA):
    Mx, sigmaA = analyticalSol(Data)
    Stress = InternalStress(Loads, WingboxArray, Data, NA)
    plt.figure(figsize=(8, 4))
    color = "tab:red"
    plt.xlabel("Span [m]")
    plt.ylabel("Lift [N]")
    plt.plot(-Data[0], Loads[3], color=color)
    format_plot()
    save_plot(".", "AppliedLoadsBeam")

    plt.show()

    plt.figure(figsize=(8, 4))
    color = "tab:red"
    plt.xlabel("Span [m]")
    plt.ylabel("Mx [Nm]")
    plt.plot((Data[0]) * -1, Loads[0], color=color, label="Model solution", linewidth=3)
    plt.plot((Data[0]) * -1, Mx, color="blue", label="Analytical solution", linewidth=3)
    plt.legend()
    format_plot()
    save_plot(".", "InternalBeam")
    plt.show()

    # # Internal stress diagram
    plt.figure(figsize=(8, 4))
    plt.plot(-Data[0], Stress[1, 0], label="Model stress distribution in point 1", linewidth=3)
    plt.plot(-Data[0], Stress[1, 1], label="Model stress distribution in point 2", linewidth=3)
    plt.plot(-Data[0], Stress[1, 2], label="Model stress distribution in point 3", linewidth=3)
    plt.plot(-Data[0], Stress[1, 3], label="Model stress distribution in point 4", linewidth=3)
    plt.plot(-Data[0], sigmaA[0], label="Analytical stress distribution in point 1", linewidth=3)
    plt.plot(-Data[0], sigmaA[1], label="Analytical stress distribution in point 2", linewidth=3)
    plt.plot(-Data[0], sigmaA[2], label="Analytical stress distribution in point 3", linewidth=3)
    plt.plot(-Data[0], sigmaA[3], label="Analytical stress distribution in point 4", linewidth=3)
    plt.legend()
    plt.ylabel("Stress [Mpa]")
    plt.xlabel("Span [m]")

    format_plot()
    save_plot(".", "StressBeam")
    plt.show()

    # Plot the wingbox structure
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axhline(NA[1], linestyle="--", color="darkgrey", zorder=1)
    ax.axvline(NA[0], linestyle="--", color="darkgrey", zorder=2)
    if indx != 2:
        ax.plot(*list(zip(*Airfoilcoordinates)), color="lightgrey", zorder=3, linewidth=3)
    ax.plot(WingboxArray[0], WingboxArray[1], color="black", zorder=4, linewidth=3)
    ax.plot(
        [WingboxArray[0, 0], WingboxArray[0, -1]],
        [WingboxArray[1, 0], WingboxArray[1, -1]],
        color="black",
        zorder=4,
        linewidth=3,
    )
    n = [1, 2, 3, 4]

    ax.scatter(WingboxArray[0], WingboxArray[1], color="orange", zorder=5)
    for i, txt in enumerate(n):
        ax.annotate(txt, (WingboxArray[0, i], WingboxArray[1, i]), zorder=6)
    plt.gca().axis("equal")
    plt.ylabel("z/MAC [-]")
    plt.xlabel("x/MAC [-]")
    format_plot()
    save_plot(".", "WingboxBeam")
    plt.show()


n = 1000
indx = 2

if indx == 0:
    Airfoilcoordinates = GetAirfoil(filenameArray[-1])
    Data = GetDataXFLR(filenameArray[indx])
    WingboxArray = GetWingbox(indx)
    NA = NeutralAxisUnit(WingboxArray)
    Loads = AppliedLoads(Data)
    Stress = InternalStress(Loads, WingboxArray, Data, NA)
    MakePlots0(indx, WingboxArray, Data, Loads, NA)
elif indx == 1:
    Airfoilcoordinates = GetAirfoil(filenameArray[-1])
    Data2 = GetDataXFLR(filenameArray[indx])
    WingboxArray = GetWingbox(indx)
    NA = NeutralAxisUnit(WingboxArray)
    Data1 = GetDataXFLR(filenameArray[indx + 1])
    Loads1 = AppliedLoads(Data1)
    Loads2 = AppliedLoads(Data2)
    Stress1 = InternalStress(Loads1, WingboxArray, Data1, NA)
    Stress2 = InternalStress(Loads2, WingboxArray, Data2, NA)
    MakePlots1(indx, WingboxArray, Data1, Data2, Loads1, Loads2, NA)
else:
    Data = ConstructDataCB(n)
    WingboxArray = GetWingbox(indx)
    NA = NeutralAxisUnit(WingboxArray)
    Loads = AppliedLoadsCB(Data)
    Loads1 = AppliedLoads(Data)

    MakePlots2(indx, WingboxArray, Data, Loads, NA)
