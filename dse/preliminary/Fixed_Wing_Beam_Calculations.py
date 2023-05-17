import matplotlib.pyplot as plt
import csv


### Cruise conditions
V_cruise = 400/3.6  # [m/s]
rho = 0.01  # [kg/m3] air density
q = 0.5 * rho * (V_cruise ** 2)  # [Pa] Dynamic pressure
n = 2.5  # [-] Maximum load factor
j = 1
if j == 0:
    ### Airfoil Data
    with open('S1223.csv') as csvfile:
        data = csv.reader(csvfile, delimiter='\t')
        Data = []
        for row in data:
            Data.append(row)
    Airfoilcoordinates = [[float(Data[i][0]), float(Data[i][1])] for i in range(len(Data))]


    x_coor_lst = [0.05, 0.5, 0.5, 0.05]  # [-] As percentage of unit chord
    y_coor_lst = [0.07, 0.12486, 0.04953, -0.01123]  # [-] As percentage of unit chord
    t_sheet_lst = [0.001, 0.01, 0.001, 0.01]  # [m] [t12, t23, t34, t41
    # b_lst = []  # [-] The length of each sheet for the unit chord: [b12, b23, b34, b41]
    # for i in range(len(x_coor_lst)-1):
    #     b_lst.append(((x_coor_lst[i]-x_coor_lst[i+1])**2+(y_coor_lst[i]-y_coor_lst[i+1])**2)**0.5)
    # b_lst.append(((x_coor_lst[0]-x_coor_lst[-1])**2+(y_coor_lst[0]-y_coor_lst[-1])**2)**0.5)
    #
    # xcg_lst = []
    # ycg_lst = []
    # for i in range(len(x_coor_lst)-1):
    #     xcg_lst.append(x_coor_lst[i] + (x_coor_lst[i + 1] - x_coor_lst[i]) / 2)
    #     ycg_lst.append(y_coor_lst[i] + (y_coor_lst[i + 1] - y_coor_lst[i]) / 2)
    # xcg_lst.append(x_coor_lst[-1] + (x_coor_lst[0] - x_coor_lst[-1]) / 2)
    # ycg_lst.append(y_coor_lst[-1] + (y_coor_lst[0] - y_coor_lst[-1]) / 2)
    #
    # def NeutralAxisUnit():
    #     xbar = sum(b_lst[i]*t_sheet_lst[i]*xcg_lst[i] for i in range(len(b_lst)))/sum(b_lst[i]*t_sheet_lst[i] for i in range(len(b_lst)))
    #     ybar = sum(b_lst[i]*t_sheet_lst[i]*ycg_lst[i] for i in range(len(b_lst)))/sum(b_lst[i]*t_sheet_lst[i] for i in range(len(b_lst)))
    #     return xbar, ybar
    #
    # def IxxCalcUnit():  # Multiply by c^4 to gain the moment of inertia at said location
    #     ybar = NeutralAxisUnit()[1]
    #     Ixx12 = (t_sheet_lst[0]*b_lst[0]**3*((y_coor_lst[1]-y_coor_lst[0]) / b_lst[0])**2)/12
    #     Ixx23 = (t_sheet_lst[1]*b_lst[1]**3)/12
    #     Ixx34 = (t_sheet_lst[2]*b_lst[2]**3*((y_coor_lst[3]-y_coor_lst[2]) / b_lst[2])**2)/12
    #     Ixx41 = (t_sheet_lst[3] * b_lst[3] ** 3) / 12
    #     IxxPart = [Ixx12, Ixx23, Ixx34, Ixx41]
    #     Ixx = sum(IxxPart[i]+t_sheet_lst[i]*b_lst[i]*(ycg_lst[i]-ybar)**2 for i in range(len(IxxPart)))
    #     return Ixx
    #
    # def IyyCalcUnit():  # Multiply by c^4 to gain the moment of inertia at said location
    #     xbar = NeutralAxisUnit()[0]
    #     Iyy12 = (t_sheet_lst[0] * b_lst[0] ** 3 * (
    #                 (x_coor_lst[1] - x_coor_lst[0]) / b_lst[0]) ** 2) / 12
    #     Iyy23 = 0
    #     Iyy34 = (t_sheet_lst[2] * b_lst[2] ** 3 * (
    #                 (x_coor_lst[3] - x_coor_lst[2]) / b_lst[2]) ** 2) / 12
    #     Iyy41 = 0
    #     IyyPart = [Iyy12, Iyy23, Iyy34, Iyy41]
    #     Iyy = sum(IyyPart[i] + t_sheet_lst[i] * b_lst[i] * (xcg_lst[i] - xbar) ** 2 for i in range(len(IyyPart)))
    #     return Iyy
    #
    # def IxyCalcUnit():  # Multiply by c^4 to gain the moment of inertia at said location
    #     xbar, ybar = NeutralAxisUnit()
    #     Ixy12 = (t_sheet_lst[0] * b_lst[0] ** 3 * (((y_coor_lst[1] - y_coor_lst[0])*(x_coor_lst[1] - x_coor_lst[0])) / (b_lst[0]**2))) / 12
    #     Ixy23 = 0
    #     Ixy34 = (t_sheet_lst[2] * b_lst[2] ** 3 * (((y_coor_lst[3] - y_coor_lst[2])*(x_coor_lst[3] - x_coor_lst[2])) / (b_lst[2]**2))) / 12
    #     Ixy41 = 0
    #     IxyPart = [Ixy12, Ixy23, Ixy34, Ixy41]
    #     Ixy = sum(IxyPart[i] + t_sheet_lst[i] * b_lst[i] * (ycg_lst[i] - ybar) * (xcg_lst[i] - xbar) for i in range(len(IxyPart)))
    #     return Ixy
else:
    # for cantilever beam
    x_coor_lst = [0.0, 1, 1, 0.0]  # [-] As percentage of unit chord
    y_coor_lst = [0.5, 0.5, -0.5, -0.5]  # [-] As percentage of unit chord
    t_sheet_lst = [0.001, 0.001, 0.001, 0.001]

b_lst = []  # [-] The length of each sheet for the unit chord: [b12, b23, b34, b41]
for i in range(len(x_coor_lst) - 1):
    b_lst.append(((x_coor_lst[i] - x_coor_lst[i + 1]) ** 2 + (y_coor_lst[i] - y_coor_lst[i + 1]) ** 2) ** 0.5)
b_lst.append(((x_coor_lst[0] - x_coor_lst[-1]) ** 2 + (y_coor_lst[0] - y_coor_lst[-1]) ** 2) ** 0.5)

xcg_lst = []
ycg_lst = []
for i in range(len(x_coor_lst) - 1):
    xcg_lst.append(x_coor_lst[i] + (x_coor_lst[i + 1] - x_coor_lst[i]) / 2)
    ycg_lst.append(y_coor_lst[i] + (y_coor_lst[i + 1] - y_coor_lst[i]) / 2)
xcg_lst.append(x_coor_lst[-1] + (x_coor_lst[0] - x_coor_lst[-1]) / 2)
ycg_lst.append(y_coor_lst[-1] + (y_coor_lst[0] - y_coor_lst[-1]) / 2)


def NeutralAxisUnit():
    xbar = sum(b_lst[i] * t_sheet_lst[i] * xcg_lst[i] for i in range(len(b_lst))) / sum(
        b_lst[i] * t_sheet_lst[i] for i in range(len(b_lst)))
    ybar = sum(b_lst[i] * t_sheet_lst[i] * ycg_lst[i] for i in range(len(b_lst))) / sum(
        b_lst[i] * t_sheet_lst[i] for i in range(len(b_lst)))
    return xbar, ybar


def IxxCalcUnit():  # Multiply by c^4 to gain the moment of inertia at said location
    ybar = NeutralAxisUnit()[1]
    Ixx12 = (t_sheet_lst[0] * b_lst[0] ** 3 * ((y_coor_lst[1] - y_coor_lst[0]) / b_lst[0]) ** 2) / 12
    Ixx23 = (t_sheet_lst[1] * b_lst[1] ** 3) / 12
    Ixx34 = (t_sheet_lst[2] * b_lst[2] ** 3 * ((y_coor_lst[3] - y_coor_lst[2]) / b_lst[2]) ** 2) / 12
    Ixx41 = (t_sheet_lst[3] * b_lst[3] ** 3) / 12
    IxxPart = [Ixx12, Ixx23, Ixx34, Ixx41]
    Ixx = sum(IxxPart[i] + t_sheet_lst[i] * b_lst[i] * (ycg_lst[i] - ybar) ** 2 for i in range(len(IxxPart)))
    return Ixx


def IyyCalcUnit():  # Multiply by c^4 to gain the moment of inertia at said location
    xbar = NeutralAxisUnit()[0]
    Iyy12 = (t_sheet_lst[0] * b_lst[0] ** 3 * (
            (x_coor_lst[1] - x_coor_lst[0]) / b_lst[0]) ** 2) / 12
    Iyy23 = 0
    Iyy34 = (t_sheet_lst[2] * b_lst[2] ** 3 * (
            (x_coor_lst[3] - x_coor_lst[2]) / b_lst[2]) ** 2) / 12
    Iyy41 = 0
    IyyPart = [Iyy12, Iyy23, Iyy34, Iyy41]
    Iyy = sum(IyyPart[i] + t_sheet_lst[i] * b_lst[i] * (xcg_lst[i] - xbar) ** 2 for i in range(len(IyyPart)))
    return Iyy


def IxyCalcUnit():  # Multiply by c^4 to gain the moment of inertia at said location
    xbar, ybar = NeutralAxisUnit()
    Ixy12 = (t_sheet_lst[0] * b_lst[0] ** 3 * (
                ((y_coor_lst[1] - y_coor_lst[0]) * (x_coor_lst[1] - x_coor_lst[0])) / (b_lst[0] ** 2))) / 12
    Ixy23 = 0
    Ixy34 = (t_sheet_lst[2] * b_lst[2] ** 3 * (
                ((y_coor_lst[3] - y_coor_lst[2]) * (x_coor_lst[3] - x_coor_lst[2])) / (b_lst[2] ** 2))) / 12
    Ixy41 = 0
    IxyPart = [Ixy12, Ixy23, Ixy34, Ixy41]
    Ixy = sum(IxyPart[i] + t_sheet_lst[i] * b_lst[i] * (ycg_lst[i] - ybar) * (xcg_lst[i] - xbar) for i in
              range(len(IxyPart)))
    return Ixy




def FlyingWing():
    ### Reading of XFLR5 Data
    with open('Flying_Wing_Data_AoA=9_3.csv') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
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

    def AppliedLoads():
        Mx_lst = [0]
        My_lst = [0]
        Mz_lst = [0]
        L_lst = [0]
        D_lst = [0]
        S_lst = []
        for i in range(len(y_lst) - 1):
            dy = y_lst[i + 1] - y_lst[i]
            c_avg = (c_lst[i + 1] + c_lst[i]) / 2
            L = q * CL_lst[i] * c_avg * dy
            D = q * (ICd_lst[i] + PCd_lst[i]) * c_avg * dy
            S_lst.append(dy * c_avg)
            L_lst.append(L)
            D_lst.append(D)
            Mx_lst.append(Mx_lst[i] + L * dy / 2 + sum(L_lst) * dy)
            My_lst.append(My_lst[i] + D * dy / 2 + sum(D_lst) * dy)
        return Mx_lst, My_lst, Mz_lst, L_lst, D_lst, S_lst

    def InternalStress():
        Mx_lst, My_lst = AppliedLoads()[0:2]
        Ixx = IxxCalcUnit()
        Iyy = IyyCalcUnit()
        Ixy = IxyCalcUnit()
        xbar, ybar = NeutralAxisUnit()
        sigmaz_1_lst = [10 ** -6 * ((Mx_lst[i] * Iyy - My_lst[i] * Ixy) * (-y_coor_lst[0] + ybar) + (
                    My_lst[i] * Ixx - Mx_lst[i] * Ixy) * (-x_coor_lst[0] + xbar)) / (
                                    (Ixx * Iyy - Ixy ** 2) * c_lst[i] ** 2) for i in range(len(y_lst))]
        sigmaz_2_lst = [10 ** -6 * ((Mx_lst[i] * Iyy - My_lst[i] * Ixy) * (-y_coor_lst[1] + ybar) + (
                    My_lst[i] * Ixx - Mx_lst[i] * Ixy) * (x_coor_lst[1] - xbar)) / (
                                    (Ixx * Iyy - Ixy ** 2) * c_lst[i] ** 2) for i in range(len(y_lst))]
        sigmaz_3_lst = [10 ** -6 *
                        ((Mx_lst[i] * Iyy - My_lst[i] * Ixy) * (-y_coor_lst[2] + ybar) + (
                                    My_lst[i] * Ixx - Mx_lst[i] * Ixy) * (-x_coor_lst[2] + xbar)) / (
                                (Ixx * Iyy - Ixy ** 2) * c_lst[i] ** 2) for i in range(len(y_lst))]
        sigmaz_4_lst = [10 ** -6 *
                        ((Mx_lst[i] * Iyy - My_lst[i] * Ixy) * (-y_coor_lst[3] + ybar) + (
                                    My_lst[i] * Ixx - Mx_lst[i] * Ixy) * (x_coor_lst[3] - xbar)) / (
                                (Ixx * Iyy - Ixy ** 2) * c_lst[i] ** 2) for i in range(len(y_lst))]
        return sigmaz_1_lst, sigmaz_2_lst, sigmaz_3_lst, sigmaz_4_lst

    ### Making the plots
    # Making the applied loads plots
    Mx_lst, My_lst, Mz_lst, L_lst, D_lst, S_lst = AppliedLoads()

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Span [m]')
    ax1.set_ylabel('Lift [N]', color=color)
    ax1.plot(y_lst, L_lst, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()

    color = 'tab:blue'
    ax2.set_ylabel('Drag [N]', color=color)  # we already handled the x-label with ax1
    ax2.plot(y_lst, D_lst, color=color, linestyle='--')
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

    # Making the internal loads plots

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Span [m]')
    ax1.set_ylabel('Mx [Nm]', color=color)
    ax1.plot(y_lst, Mx_lst, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()

    color = 'tab:blue'
    ax2.set_ylabel('My [Nm]', color=color)  # we already handled the x-label with ax1
    ax2.plot(y_lst, My_lst, color=color, linestyle='--')
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

    # Internal stress diagram
    sigmaz_1_lst, sigmaz_2_lst, sigmaz_3_lst, sigmaz_4_lst = InternalStress()
    plt.plot(y_lst, sigmaz_1_lst, label='Stress distribution in point 1')
    plt.plot(y_lst, sigmaz_2_lst, label='Stress distribution in point 2')
    plt.plot(y_lst, sigmaz_3_lst, label='Stress distribution in point 3')
    plt.plot(y_lst, sigmaz_4_lst, label='Stress distribution in point 4')
    plt.legend()
    plt.ylabel('Stress [Mpa]')
    plt.xlabel('Span [m]')
    plt.show()

    # Plot the wingbox structure
    fig, ax = plt.subplots()
    xbar, ybar = NeutralAxisUnit()

    ax.axhline(ybar, linestyle='--', color='darkgrey', zorder=1)
    ax.axvline(xbar, linestyle='--', color='darkgrey', zorder=2)
    ax.plot(*list(zip(*Airfoilcoordinates)), color='lightgrey', zorder=3)
    ax.plot(x_coor_lst, y_coor_lst, color='black', zorder=4)
    ax.plot([x_coor_lst[0], x_coor_lst[-1]], [y_coor_lst[0], y_coor_lst[-1]], color='black', zorder=4)
    n = [1, 2, 3, 4]

    ax.scatter(x_coor_lst, y_coor_lst, color='orange', zorder=5)
    for i, txt in enumerate(n):
        ax.annotate(txt, (x_coor_lst[i], y_coor_lst[i]), zorder=6)
    plt.gca().axis("equal")
    plt.ylabel('z/MAC [-]')
    plt.xlabel('x/MAC [-]')
    plt.show()

def BiWing():
    ### Reading of XFLR5 Data
    y_lst = []
    c_lst = []
    CL_lst = []
    PCd_lst = []
    ICd_lst = []
    with open('Bi_Wing_1_a=14.00_v=111.00ms.csv') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        Data = []
        for row in data:
            Data.append(row)

    y_lst.append([float(i[0]) for i in Data[1:]]) # [m] Y-step of the applied loads b
    c_lst.append([float(i[1]) for i in Data[1:]])  # [m] Y-step of the applied loads b
    AoA_lst = [float(i[2]) for i in Data[1:]]  # [degrees] induced angle of angle of attack
    CL_lst.append([float(i[3]) for i in Data[1:]])  # [-] Lift coefficient of the wing
    PCd_lst.append([float(i[4]) for i in Data[1:]])  # [-] Parasite drag coefficient of the wing
    ICd_lst.append([float(i[5]) for i in Data[1:]])  # [-] Induced drag coefficient of the wing
    CmGeom_lst = [float(i[6]) for i in Data[1:]]  # [-] Moment coefficient
    CmAirfchord4_lst = [float(i[7]) for i in Data[1:]]  # [-] Moment coefficient at quarter chord

    with open('Bi_Wing_2_a=14.00_v=111.00ms.csv') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        Data = []
        for row in data:
            Data.append(row)

    y_lst.append([float(i[0]) for i in Data[1:]])  # [m] Y-step of the applied loads b
    c_lst.append([float(i[1]) for i in Data[1:]])  # [m] Y-step of the applied loads b
    AoA_lst = [float(i[2]) for i in Data[1:]]  # [degrees] induced angle of angle of attack
    CL_lst.append([float(i[3]) for i in Data[1:]])  # [-] Lift coefficient of the wing
    PCd_lst.append([float(i[4]) for i in Data[1:]])  # [-] Parasite drag coefficient of the wing
    ICd_lst.append([float(i[5]) for i in Data[1:]])  # [-] Induced drag coefficient of the wing
    CmGeom_lst = [float(i[6]) for i in Data[1:]]  # [-] Moment coefficient
    CmAirfchord4_lst = [float(i[7]) for i in Data[1:]]  # [-] Moment coefficient at quarter chord

    def AppliedLoads():
        Mx_lst = []
        My_lst = []
        Mz_lst = []
        L_lst = []
        D_lst = []
        S_lst = []
        for j in range(len(y_lst)):
            Mxi_lst = [0]
            Myi_lst = [0]
            Mzi_lst = [0]
            Li_lst = [0]
            Di_lst = [0]
            Si_lst = []
            for i in range(len(y_lst[j]) - 1):
                dy = y_lst[j][i + 1] - y_lst[j][i]
                c_avg = (c_lst[j][i + 1] + c_lst[j][i]) / 2
                L = q * CL_lst[j][i] * c_avg * dy
                D = q * (ICd_lst[j][i] + PCd_lst[j][i]) * c_avg * dy
                Si_lst.append(dy * c_avg)
                Li_lst.append(L)
                Di_lst.append(D)
                Mxi_lst.append(Mxi_lst[i] + L * dy / 2 + sum(Li_lst) * dy)
                Myi_lst.append(Myi_lst[i] + D * dy / 2 + sum(Di_lst) * dy)
            L_lst.append(Li_lst)
            D_lst.append(Di_lst)
            Mx_lst.append(Mxi_lst)
            My_lst.append(Myi_lst)
        return Mx_lst, My_lst, Mz_lst, L_lst, D_lst, S_lst

    def InternalStress():
        Mx_lst, My_lst = AppliedLoads()[0:2]
        Ixx = IxxCalcUnit()
        Iyy = IyyCalcUnit()
        Ixy = IxyCalcUnit()
        xbar, ybar = NeutralAxisUnit()
        sigmaz_1_lst = []
        sigmaz_2_lst = []
        sigmaz_3_lst =[]
        sigmaz_4_lst = []
        for j in range(len(y_lst)):
            sigmazi_1_lst = [10 ** -6 * ((Mx_lst[j][i] * Iyy - My_lst[j][i] * Ixy) * (-y_coor_lst[0] + ybar) + (
                    My_lst[j][i] * Ixx - Mx_lst[j][i] * Ixy) * (-x_coor_lst[0] + xbar)) / (
                                    (Ixx * Iyy - Ixy ** 2) * c_lst[j][i] ** 2) for i in range(len(y_lst[j]))]
            sigmazi_2_lst = [10 ** -6 * ((Mx_lst[j][i] * Iyy - My_lst[j][i] * Ixy) * (-y_coor_lst[1] + ybar) + (
                    My_lst[j][i] * Ixx - Mx_lst[j][i] * Ixy) * (x_coor_lst[1] - xbar)) / (
                                    (Ixx * Iyy - Ixy ** 2) * c_lst[j][i] ** 2) for i in range(len(y_lst[j]))]
            sigmazi_3_lst = [10 ** -6 *
                            ((Mx_lst[j][i] * Iyy - My_lst[j][i] * Ixy) * (-y_coor_lst[2] + ybar) + (
                                    My_lst[j][i] * Ixx - Mx_lst[j][i] * Ixy) * (-x_coor_lst[2] + xbar)) / (
                                    (Ixx * Iyy - Ixy ** 2) * c_lst[j][i] ** 2) for i in range(len(y_lst[j]))]
            sigmazi_4_lst = [10 ** -6 *
                            ((Mx_lst[j][i] * Iyy - My_lst[j][i] * Ixy) * (-y_coor_lst[3] + ybar) + (
                                    My_lst[j][i] * Ixx - Mx_lst[j][i] * Ixy) * (x_coor_lst[3] - xbar)) / (
                                    (Ixx * Iyy - Ixy ** 2) * c_lst[j][i] ** 2) for i in range(len(y_lst[j]))]
            sigmaz_1_lst.append(sigmazi_1_lst)
            sigmaz_2_lst.append(sigmazi_2_lst)
            sigmaz_3_lst.append(sigmazi_3_lst)
            sigmaz_4_lst.append(sigmazi_4_lst)
        return sigmaz_1_lst, sigmaz_2_lst, sigmaz_3_lst, sigmaz_4_lst

    ### Making the plots
    # Making the applied loads plots
    Mx_lst, My_lst, Mz_lst, L_lst, D_lst, S_lst = AppliedLoads()


    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Span [m]')
    ax1.set_ylabel('Lift [N]', color=color)
    ax1.plot(y_lst[0], L_lst[0], color=color,linestyle='--', label='Lift lower wing')
    ax1.plot(y_lst[1], L_lst[1], color=color, label='Lift upper wing')
    ax1.tick_params(axis='y', labelcolor=color)
    plt.legend()
    ax2 = ax1.twinx()

    color = 'tab:blue'
    ax2.set_ylabel('Drag [N]', color=color)  # we already handled the x-label with ax1
    ax2.plot(y_lst[0], D_lst[0], color=color, linestyle='--', label='Drag lower wing')
    ax2.plot(y_lst[1], D_lst[1], color=color, label='Drag upper wing')
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.legend()
    plt.show()

    # Making the internal loads plots

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Span [m]')
    ax1.set_ylabel('Mx [Nm]', color=color)
    ax1.plot(y_lst[0], Mx_lst[0], color=color, linestyle='--', label='Lift lower wing')
    ax1.plot(y_lst[1], Mx_lst[1], color=color, label='Lift upper wing')
    ax1.tick_params(axis='y', labelcolor=color)
    plt.legend()
    ax2 = ax1.twinx()

    color = 'tab:blue'
    ax2.set_ylabel('My [Nm]', color=color)  # we already handled the x-label with ax1
    ax2.plot(y_lst[0], My_lst[0], color=color, linestyle='--', label='Drag lower wing')
    ax2.plot(y_lst[1], My_lst[1], color=color, label='Drag upper wing')
    ax2.tick_params(axis='y', labelcolor=color)


    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

    # Internal stress diagram
    sigmaz_1_lst, sigmaz_2_lst, sigmaz_3_lst, sigmaz_4_lst = InternalStress()
    print(sigmaz_1_lst)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(y_lst[0], sigmaz_1_lst[0], label='Stress distribution in point 1')
    ax1.plot(y_lst[0], sigmaz_2_lst[0], label='Stress distribution in point 2')
    ax1.plot(y_lst[0], sigmaz_3_lst[0], label='Stress distribution in point 3')
    ax1.plot(y_lst[0], sigmaz_4_lst[0], label='Stress distribution in point 4')
    ax1.legend()
    plt.ylabel('Stress [Mpa]')
    plt.xlabel('Span [m]')
    ax2.plot(y_lst[1], sigmaz_1_lst[1], label='Stress distribution in point 1')
    ax2.plot(y_lst[1], sigmaz_2_lst[1], label='Stress distribution in point 2')
    ax2.plot(y_lst[1], sigmaz_3_lst[1], label='Stress distribution in point 3')
    ax2.plot(y_lst[1], sigmaz_4_lst[1], label='Stress distribution in point 4')
    ax2.legend()
    plt.ylabel('Stress [Mpa]')
    plt.xlabel('Span [m]')
    plt.show()


    # Plot the wingbox structure
    fig, ax = plt.subplots()
    xbar, ybar = NeutralAxisUnit()

    ax.axhline(ybar, linestyle='--', color='darkgrey', zorder=1)
    ax.axvline(xbar, linestyle='--', color='darkgrey', zorder=2)
    ax.plot(*list(zip(*Airfoilcoordinates)), color='lightgrey', zorder=3)
    ax.plot(x_coor_lst, y_coor_lst, color='black', zorder=4)
    ax.plot([x_coor_lst[0], x_coor_lst[-1]], [y_coor_lst[0], y_coor_lst[-1]], color='black', zorder=4)
    n = [1, 2, 3, 4]

    ax.scatter(x_coor_lst, y_coor_lst, color='orange', zorder=5)
    for i, txt in enumerate(n):
        ax.annotate(txt, (x_coor_lst[i], y_coor_lst[i]), zorder=6)
    plt.gca().axis("equal")
    plt.ylabel('z/MAC [-]')
    plt.xlabel('x/MAC [-]')
    plt.show()

def CantileverBeam():
    ### Reading of XFLR5 Data
    with open('Cantilever_Beam_More_Points.csv') as csvfile:
        data = csv.reader(csvfile, delimiter=';')
        Data = []
        for row in data:
            Data.append(row)

    y_lst = [float(i[0]) for i in Data[1:]]  # [m] Y-step of the applied loads b
    c_lst = [float(i[1]) for i in Data[1:]]  # [m] Y-step of the applied loads b
    AoA_lst = [float(i[2]) for i in Data[1:]]  # [degrees] induced angle of angle of attack
    CL_lst = [float(i[3]) for i in Data[1:]]  # [-] Lift coefficient of the wing
    PCd_lst = [float(i[4]) for i in Data[1:]]  # [-] Parasite drag coefficient of the wing
    ICd_lst = [float(i[5]) for i in Data[1:]]  # [-] Induced drag coefficient of the wing

    def AppliedLoads():
        Mx_lst = [0]
        My_lst = [0]
        Mz_lst = [0]
        L_lst = [0]
        D_lst = [0]
        S_lst = []
        w_lst = []
        for i in range(len(y_lst) - 1):
            dy = y_lst[i + 1] - y_lst[i]
            c_avg = (c_lst[i + 1] + c_lst[i]) / 2
            L = 300
            w = 300*dy
            w_lst.append(w)
            D = q * (ICd_lst[i] + PCd_lst[i]) * c_avg * dy
            S_lst.append(dy * c_avg)
            L_lst.append(L)
            D_lst.append(D)
            Mx_lst.append(Mx_lst[i] + w * dy / 2+sum(w_lst)*dy)
            My_lst.append(My_lst[i] + D * dy / 2)
        print(Mx_lst[-1])
        return Mx_lst, My_lst, Mz_lst, L_lst, D_lst, S_lst

    def InternalStress():
        Mx_lst, My_lst = AppliedLoads()[0:2]
        Ixx = IxxCalcUnit()
        Iyy = IyyCalcUnit()
        Ixy = IxyCalcUnit()
        print(Ixy)
        xbar, ybar = NeutralAxisUnit()

        sigmaz_1_lst = [10 ** -6 * ((Mx_lst[i] * Iyy - My_lst[i] * Ixy) * (-y_coor_lst[0] + ybar) + (
                My_lst[i] * Ixx - Mx_lst[i] * Ixy) * (-x_coor_lst[0] + xbar)) / (
                                (Ixx * Iyy - Ixy ** 2) * c_lst[i] ** 2) for i in range(len(y_lst))]
        sigmaz_2_lst = [10 ** -6 * ((Mx_lst[i] * Iyy - My_lst[i] * Ixy) * (-y_coor_lst[1] + ybar) + (
                My_lst[i] * Ixx - Mx_lst[i] * Ixy) * (x_coor_lst[1] - xbar)) / (
                                (Ixx * Iyy - Ixy ** 2) * c_lst[i] ** 2) for i in range(len(y_lst))]
        sigmaz_3_lst = [10 ** -6 *
                        ((Mx_lst[i] * Iyy - My_lst[i] * Ixy) * (-y_coor_lst[2] + ybar) + (
                                My_lst[i] * Ixx - Mx_lst[i] * Ixy) * (-x_coor_lst[2] + xbar)) / (
                                (Ixx * Iyy - Ixy ** 2) * c_lst[i] ** 2) for i in range(len(y_lst))]
        sigmaz_4_lst = [10 ** -6 *
                        ((Mx_lst[i] * Iyy - My_lst[i] * Ixy) * (-y_coor_lst[3] + ybar) + (
                                My_lst[i] * Ixx - Mx_lst[i] * Ixy) * (x_coor_lst[3] - xbar)) / (
                                (Ixx * Iyy - Ixy ** 2) * c_lst[i] ** 2) for i in range(len(y_lst))]
        print(max(sigmaz_4_lst))
        return sigmaz_1_lst, sigmaz_2_lst, sigmaz_3_lst, sigmaz_4_lst



    ### Making the plots
    # Making the applied loads plots
    Mx_lst, My_lst, Mz_lst, L_lst, D_lst, S_lst = AppliedLoads()

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Span [m]')
    ax1.set_ylabel('Lift [N]', color=color)
    ax1.plot(y_lst, L_lst, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()

    # color = 'tab:blue'
    # ax2.set_ylabel('Drag [N]', color=color)  # we already handled the x-label with ax1
    # ax2.plot(y_lst, D_lst, color=color, linestyle='--')
    # ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

    # Making the internal loads plots

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Span [m]')
    ax1.set_ylabel('Mx [Nm]', color=color)
    ax1.plot(y_lst, Mx_lst, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()

    # color = 'tab:blue'
    # ax2.set_ylabel('My [Nm]', color=color)  # we already handled the x-label with ax1
    # ax2.plot(y_lst, My_lst, color=color, linestyle='--')
    # ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

    # Internal stress diagram
    sigmaz_1_lst, sigmaz_2_lst, sigmaz_3_lst, sigmaz_4_lst = InternalStress()
    plt.plot(y_lst, sigmaz_1_lst, label='Stress distribution in point 1')
    plt.plot(y_lst, sigmaz_2_lst, label='Stress distribution in point 2')
    plt.plot(y_lst, sigmaz_3_lst, label='Stress distribution in point 3')
    plt.plot(y_lst, sigmaz_4_lst, label='Stress distribution in point 4')
    plt.legend()
    plt.ylabel('Stress [Mpa]')
    plt.xlabel('Span [m]')
    plt.show()

    # Plot the wingbox structure
    fig, ax = plt.subplots()
    xbar, ybar = NeutralAxisUnit()

    ax.axhline(ybar, linestyle='--', color='darkgrey', zorder=1)
    ax.axvline(xbar, linestyle='--', color='darkgrey', zorder=2)
    ax.plot(x_coor_lst, y_coor_lst, color='black', zorder=4)
    ax.plot([x_coor_lst[0], x_coor_lst[-1]], [y_coor_lst[0], y_coor_lst[-1]], color='black', zorder=4)
    n = [1, 2, 3, 4]

    ax.scatter(x_coor_lst, y_coor_lst, color='orange', zorder=5)
    for i, txt in enumerate(n):
        ax.annotate(txt, (x_coor_lst[i], y_coor_lst[i]), zorder=6)
    plt.gca().axis("equal")
    plt.ylabel('z/MAC [-]')
    plt.xlabel('x/MAC [-]')
    plt.show()



# 0 for flying wing, 1 for conventional fixed wing aircraft, 2 for Cantilever Beam
i = 2
if i == 0:
    FlyingWing()
elif i == 1:
    BiWing()
else:
    CantileverBeam()



