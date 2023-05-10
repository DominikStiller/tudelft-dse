import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipkinc, ellipeinc

## Constants
rho0 = 1.0 * 10 ** -2  # [kg/m3] Zero elevation density of Mars
T0_low = 145  # [K] Lowest zero elevation temperature
T0 = 273 # [K] T0 per definition
T0_high = 284  # [K] Highest zero elevation temperature
S = 222  # [K] Entropy
mu0 = 1.37*10**-5  # [Ns/m2]
R_M = 191  # [J/K Kg] Specific Gas constant of Martian Air
R_universal = 8.314  # [J/mol] Universal gas constant
y = 1.2  # [-] Ratio of specific heats
g_M = 3.71  # [m/s2] Gravitational acceleration of Mars
M_H2 = 2.02 * 10 ** -3  # [kg/mol] Molar mass of hydrogen
sigma_yield_kevlar = 1240 * 10 ** 6  # Yield strength kevlar
sigma_yield_polyethylene = 38 * 10 ** 6 # Yield strength polyethylene
p_ratio = 1.005

R_H2 = R_universal/M_H2

p0_low = rho0*T0_low*R_M
rho_H2_low = p_ratio*p0_low/(R_H2*T0_low)
drho_low = rho0-rho_H2_low

p0_high = rho0*T0_high*R_M
rho_H2_high = p_ratio*p0_high/(R_H2*T0_high)
drho_high = rho0-rho_H2_high

drho = drho_low  # Choose lowest drho value
p0 = p0_low

m_goal = 2700
t_min = 0.000101
m_lst = []
j_lst = []

def mass_convergence(m_misc):
    m = m_misc
    for j in range(20):
        j_lst.append(j)
        V = m/drho
        r = ((3)/(4*np.pi)*V)**(1/3)
        t_req = (p_ratio - 1) * p0 * r / (2 * sigma_yield_polyethylene)
        if t_req <= t_min:
            t = t_min
        else:
            t = t_req
        S = 4 * np.pi * r**2
        rho_kevlar = 1380 # [kg/m3]
        rho_polyethylene = 940  # [kg/m^3]
        m = S*t*rho_polyethylene + m_misc
        m_lst.append(m)
    return m, m_lst, j_lst, V, r, S

for m_misc in list(np.arange(1,3010,1)):
    m_end, m_lst, j_lst, V, r, S = mass_convergence(m_misc)
    if m_end >= m_goal:
        print(f"\nFinal total mass = {m_end} kg")
        print(f"Effective payload = {m_misc} kg\n")
        print(f"Total volume = {V} m^3 and radius of sphere = {r} m with surface area = {S} m^2")
        break

# plt.plot(j_lst, m_lst)
# plt.show()

def Reynolds():
    D_m = ((3 * 8 * V) / (4 * np.pi * F_ratio)) ** (1 / 3)
    L_a = F_ratio * D_m
    mu = mu0*(T0_low/T0)**(3/2)*((T0 + S)/(T0_low+S))
    Re = (rho0*V_cruise*L_a)/(mu)
    return Re, L_a, D_m

# CF = 0.08  # Skin friction coefficient of Ultra-High-Molecular-Weight Polyethylene
F_ratio = 4.65  # Fineness ration L_a / D_m
V_cruise = 400/3.6

Re, L_a, D_m = Reynolds()
CF = 0.045*(Re**(-1/6))  # Skin friction coefficient by literature

# a = c = D_m/2
# b = L_a/2
#
# phi = np.arccos(c/a)
# m = (a**2 * (b**2 - c**2)) / (b**2 * (a**2 - c**2))
#
# temp = ellipeinc(phi, m)*np.sin(phi)**2 + ellipkinc(phi, m)*np.cos(phi)**2
# ellipsoid_area = 2*np.pi*(c**2 + a*b*temp/np.sin(phi))
# print(f"Ellipsoid surface area = {ellipsoid_area}")

CD0 = CF * (4*(F_ratio)**(1/3) + 6*(F_ratio)**(-1.2) + 24*(F_ratio)**(-2.7))
print(f"New hull length = {L_a} m")
print(f"New hull diameter = {D_m} m")
print(f"Design Reynolds number = {Re}")
print(f"Design skin friction coefficient = {CF}")
print(f"Design CD0 = {CD0}")

CD = CD0
D = 0.5*rho0*(V_cruise**2)*np.pi*((D_m/2)**2)*CD
print(f"Drag force = {D} N")

### Hybrid Blimp-Wing Design ###
m_total = m_end
# m_total = m_wing + m_balloonskin + m_misc2

L_wing = m_total*g_M - V*drho*g_M

def S_wing(L_wing):
    CL = 1.5
    S = L_wing/(1/2 * rho0 * V_cruise ** 2 * CL)
    return S

def m_wing(S):
    wingloading = 848.8/179
    m_wing = S * wingloading
    return m_wing

def mass_convergence2(m_wing,m_misc):
    m = m_wing + m_misc
    for j in range(20):
        j_lst.append(j)
        V = m/drho
        r = ((3)/(4*np.pi)*V)**(1/3)
        t_req = (p_ratio - 1) * p0 * r / (2 * sigma_yield_polyethylene)
        if t_req <= t_min:
            t = t_min
        else:
            t = t_req
        S = 4 * np.pi * r**2
        rho_kevlar = 1380 # [kg/m3]
        rho_polyethylene = 940  # [kg/m^3]
        m += S*t*rho_polyethylene
        m_lst.append(m)
    return m, m_lst, j_lst, V, r, S

for m_misc in list(np.arange(1,3010,1)):
    m_end, m_lst, j_lst, V, r, S = mass_convergence(m_misc)
    if m_end >= m_goal:
        print(f"\nFinal total mass = {m_end} kg")
        print(f"Effective payload = {m_misc} kg\n")
        print(f"Total volume = {V} m^3 and radius of sphere = {r} m with surface area = {S} m^2")
        break

CD = CD0_blimp + CD0_wing + CL**2/(np.pi*A*e)