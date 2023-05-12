import matplotlib.pyplot as plt
import csv


### Reading of XFLR5 Data
with open('Fying_Wing_Data_AoA=9.csv') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    Data = []
    for row in data:
        Data.append(row)

y_lst = [float(i[0]) for i in Data[1:]]  # [m] Y-step of the applied loads b
AoA_lst = [float(i[1]) for i in Data[1:]]  # [degrees] induced angle of angle of attack
CL_lst = [float(i[2]) for i in Data[1:]]  # [-] Lift coefficient of the wing
PCd_lst = [float(i[3]) for i in Data[1:]]  # [-] Parasite drag coefficient of the wing
ICd_lst = [float(i[4]) for i in Data[1:]]  # [-] Induced drag coefficient of the wing
CmGeom_lst = [float(i[5]) for i in Data[1:]]  # [-] Moment coefficient
CmAirfchord4_lst = [float(i[6]) for i in Data[1:]]  # [-] Moment coefficient at quarter chord

### Wing Planform




x_coor_lst = [0.05, 0.5, 0.5, 0.05]  # [-] As percentage of unit chord
y_coor_lst = [0.02292, 0.12486, 0.04953, -0.01123]  # [-] As percentage of unit chord
t_sheet_lst = [0.001, 0.01, 0.001, 0.01]  # [m] [t12, t23, t34, t41
b_lst = []  # [-] The length of each sheet for the unit chord: [b12, b23, b34, b41]
for i in range(len(x_coor_lst)-1):
    b_lst.append(((x_coor_lst[i]-x_coor_lst[i+1])**2+(y_coor_lst[i]-y_coor_lst[i+1])**2)**0.5)
b_lst.append(((x_coor_lst[0]-x_coor_lst[-1])**2+(y_coor_lst[0]-y_coor_lst[-1])**2)**0.5)

xcg_lst = []
ycg_lst = []
for i in range(len(x_coor_lst)-1):
    xcg_lst.append(x_coor_lst[i] + (x_coor_lst[i + 1] - x_coor_lst[i]) / 2)
    ycg_lst.append(y_coor_lst[i] + (y_coor_lst[i + 1] - y_coor_lst[i]) / 2)
xcg_lst.append(x_coor_lst[-1] + (x_coor_lst[-1] - x_coor_lst[0]) / 2)
ycg_lst.append(y_coor_lst[-1] + (y_coor_lst[-1] - y_coor_lst[0]) / 2)

with open('Fying_Wing_Data_AoA=9.csv') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    Data = []
    for row in data:
        Data.append(row)

y_lst = [float(i[0]) for i in Data[1:]]  # [m] Y-step of the applied loads b
AoA_lst = [float(i[1]) for i in Data[1:]]  # [degrees] induced angle of angle of attack
CL_lst = [float(i[2]) for i in Data[1:]]  # [-] Lift coefficient of the wing
PCd_lst = [float(i[3]) for i in Data[1:]]  # [-] Parasite drag coefficient of the wing
ICd_lst = [float(i[4]) for i in Data[1:]]  # [-] Induced drag coefficient of the wing
CmGeom_lst = [float(i[5]) for i in Data[1:]]  # [-] Moment coefficient
CmAirfchord4_lst = [float(i[6]) for i in Data[1:]]  # [-] Moment coefficient at quarter chord


#### Wing geometry Flying Wing
cr = 5.6167  # [m]
ct = 1.123345  # [m]
b = 67.4  # [m]


def chord(y):
    c = cr - (((cr-ct)*2)/b) * y
    return c

def AppliedLoads():
    with open('Fying_Wing_Data_AoA=9.csv') as csvfile:
        data = csv.reader(csvfile, delimiter=';')
        Data = []
        for row in data:
            Data.append(row)
    print
    return Mx_lst, My_lst, Mz_lst




def NeutralAxisUnit():
    xbar = sum(b_lst[i]*t_sheet_lst[i]*xcg_lst[i] for i in range(len(b_lst)))/sum(b_lst[i]*t_sheet_lst[i] for i in range(len(b_lst)))
    ybar = sum(b_lst[i]*t_sheet_lst[i]*ycg_lst[i] for i in range(len(b_lst)))/sum(b_lst[i]*t_sheet_lst[i] for i in range(len(b_lst)))
    return xbar, ybar


def IxxCalcUnit():  # Multiply by c^4 to gain the moment of inertia at said location
    ybar = NeutralAxisUnit()[1]
    Ixx12 = (t_sheet_lst[0]*b_lst[0]**3*((y_coor_lst[1]-y_coor_lst[0]) / b_lst[0])**2)/12
    Ixx23 = (t_sheet_lst[1]*b_lst[1]**3)/12
    Ixx34 = (t_sheet_lst[2]*b_lst[2]**3*((y_coor_lst[3]-y_coor_lst[2]) / b_lst[2])**2)/12
    Ixx41 = (t_sheet_lst[3] * b_lst[3] ** 3) / 12
    IxxPart = [Ixx12, Ixx23, Ixx34, Ixx41]
    Ixx = sum(IxxPart[i]+t_sheet_lst[i]*b_lst[i]*(ycg_lst[i]-ybar)**2 for i in range(len(IxxPart)))
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
    Ixy12 = (t_sheet_lst[0] * b_lst[0] ** 3 * (((y_coor_lst[1] - y_coor_lst[0])*(x_coor_lst[1] - x_coor_lst[0])) / (b_lst[0]**2))) / 12
    Ixy23 = 0
    Ixy34 = (t_sheet_lst[2] * b_lst[2] ** 3 * (((y_coor_lst[3] - y_coor_lst[2])*(x_coor_lst[3] - x_coor_lst[2])) / (b_lst[2]**2))) / 12
    Ixy41 = 0
    IxyPart = [Ixy12, Ixy23, Ixy34, Ixy41]
    Ixy = sum(IxyPart[i] + t_sheet_lst[i] * b_lst[i] * (ycg_lst[i] - ybar) * (xcg_lst[i] - xbar) for i in range(len(IxyPart)))
    return Ixy

def InternalStress():
    Mx_lst, My_lst = AppliedLoads()[0:1]
    Ixx = IxxCalcUnit()
    Iyy = IyyCalcUnit()
    Ixy = IxyCalcUnit()

    sigmaz_1_lst = [((Mx_lst[i]*Iyy-My_lst[i]*Ixy)*y_coor_lst[0]+(My_lst[i]*Ixx-Mx_lst*Ixy)*x_coor_lst[0])/((Ixx*Iyy-Ixy**2)*chord(y_lst[i])**2) for i in range(len(y_lst))]
    sigmaz_2_lst = [((Mx_lst[i]*Iyy-My_lst[i]*Ixy)*y_coor_lst[1]+(My_lst[i]*Ixx-Mx_lst*Ixy)*x_coor_lst[1])/((Ixx*Iyy-Ixy**2)*chord(y_lst[i])**2) for i in range(len(y_lst))]
    sigmaz_3_lst = [
        ((Mx_lst[i] * Iyy - My_lst[i] * Ixy) * y_coor_lst[2] + (My_lst[i] * Ixx - Mx_lst * Ixy) * x_coor_lst[2]) / (
                    (Ixx * Iyy - Ixy ** 2) * chord(y_lst[i]) ** 2) for i in range(len(y_lst))]
    sigmaz_4_lst = [
        ((Mx_lst[i] * Iyy - My_lst[i] * Ixy) * y_coor_lst[3] + (My_lst[i] * Ixx - Mx_lst * Ixy) * x_coor_lst[3]) / (
                    (Ixx * Iyy - Ixy ** 2) * chord(y_lst[i]) ** 2) for i in range(len(y_lst))]
    return sigmaz_1_lst, sigmaz_2_lst, sigmaz_3_lst, sigmaz_4_lst




plt.plot(x_coor_lst, y_coor_lst)
plt.plot([x_coor_lst[0], x_coor_lst[-1]], [y_coor_lst[0], y_coor_lst[-1]])
plt.show()