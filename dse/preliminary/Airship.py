import numpy as np
import scipy as sp

## Constants
rho0 = 1.5 * 10 ** -2  # [kg/m3] Zero elevation density of Mars
T0_low = 145  # [K] Lowest zero elevation temperature
T0_high = 284  # [K] Highest zero elevation temperature
T0_Earth = 273  # [K] T0 per definition
S = 222  # [K] Entropy
mu0 = 1.37*10**-5  # [Ns/m2]
R_M = 191  # [J/K Kg] Specific Gas constant of Martian Air
R_universal = 8.314  # [J/mol] Universal gas constant
y_M = 1.2  # [-] Ratio of specific heats on Mars
g_M = 3.71  # [m/s2] Gravitational acceleration of Mars
M_H2 = 2.02 * 10 ** -3  # [kg/mol] Molar mass of hydrogen
sigma_yield_kevlar = 1240 * 10 ** 6  # Yield strength kevlar
sigma_yield_polyethylene = 38 * 10 ** 6  # Yield strength polyethylene
p_ratio = 1.005  # Pressurization ratio between outside and inside pressure

R_H2 = R_universal/M_H2  # Specific gas constant of hydrogen on Mars

# Density difference in low temperature atmosphere
p0_low = rho0*T0_low*R_M
rho_H2_low = p_ratio*p0_low/(R_H2*T0_low)
drho_low = rho0-rho_H2_low

# Density difference in high temperature atmosphere
p0_high = rho0*T0_high*R_M
rho_H2_high = p_ratio*p0_high/(R_H2*T0_high)
drho_high = rho0-rho_H2_high

drho = drho_low  # Choose smallest difference in density as most constraining
p0 = p0_low

m_design = 2700  # [kg] Mass to be designed for
t_min = 0.0001  # [m] Smallest achievable skin thickness

# Material properties
rho_kevlar = 1380  # [kg/m^3]
rho_polyethylene = 940  # [kg/m^3]
sigma_yield_kevlar = 1240 * 10 ** 6  # [Pa]
sigma_yield_polyethylene = 38 * 10 ** 6  # [Pa]

# Give in effective payload mass
# If new mass is above design mass, effective payload mass is too high
def mass_update(m_misc):
    m_skin = m_design - m_misc
    S_skin = (m_skin) / (rho_polyethylene * t_min)
    r = np.sqrt((S_skin / (4 * np.pi)))
    V_new = (4 / 3) * np.pi * r ** 3
    m_new = drho * V_new
    return m_new, V_new, S_skin, r

def residual(m_misc):
    return mass_update(m_misc)[0] - m_design

res = sp.optimize.least_squares(residual, x0=10)
m_misc_final = res.x[0]
print(f"#########")
print(f"Blimp only optimization")
print(f"#########\n")
print(f"Heighest effective mass without exceeding MTOM =  {m_misc_final} [kg]")
print(f"Total MTOM = {mass_update(m_misc_final)[0]} [kg]")
print(f"Volume is then {mass_update(m_misc_final)[1]} [m^3]")
print(f"Skin surface area is then {mass_update(m_misc_final)[2]} [m^2]")
print(f"Radius of blimp is then {mass_update(m_misc_final)[3]} [m]")

def Reynolds(V,F_ratio,V_cruise):
    D_m = ((3 * 8 * V) / (4 * np.pi * F_ratio)) ** (1 / 3)
    L_a = F_ratio * D_m
    mu = mu0 * (T0_low / T0_Earth) ** (3 / 2) * ((T0_Earth + S) / (T0_low + S))
    Re = (rho0 * V_cruise * L_a) / (mu)
    return Re, L_a, D_m

F_ratio = 4.65  # Fineness ration L_a / D_m
V_cruise = 100 / 3.6

Re, L_a, D_m = Reynolds(mass_update(m_misc_final)[1],F_ratio,V_cruise)

CF = 0.045 * (Re ** (-1 / 6))  # Skin friction coefficient by literature
CD0 = CF * (4 * (F_ratio) ** (1 / 3) + 6 * (F_ratio) ** (-1.2) + 24 * (F_ratio) ** (-2.7))
CD = CD0
D = 0.5 * rho0 * (V_cruise ** 2) * np.pi * ((D_m / 2) ** 2) * CD

print(f"New hull length = {L_a} [m]")
print(f"New hull diameter = {D_m} [m]")
print(f"Design Reynolds number = {Re}")
print(f"Design skin friction coefficient = {CF}")
print(f"Design CD0 = {CD0}")
print(f"Drag force = {D} [N]\n")

print(f"#########")
print(f"Blimp drag optimization")
print(f"#########")

def blimp_drag(m_misc):
    m_skin = m_design - m_misc
    S_skin = (m_skin) / (rho_polyethylene * t_min)
    r = np.sqrt((S_skin / (4 * np.pi)))
    V_new = (4 / 3) * np.pi * r ** 3
    m_new = drho * V_new
    Re, L_a, D_m = Reynolds(V_new, F_ratio, V_cruise)
    CF = 0.045 * (Re ** (-1 / 6))  # Skin friction coefficient by literature
    CD0 = CF * (4 * (F_ratio) ** (1 / 3) + 6 * (F_ratio) ** (-1.2) + 24 * (F_ratio) ** (-2.7))
    CD = CD0
    D = 0.5 * rho0 * (V_cruise ** 2) * np.pi * ((D_m / 2) ** 2) * CD
    return D

res = sp.optimize.least_squares(blimp_drag, x0=10, bounds=(-np.inf, m_misc_final))
print(f"Lowest drag found at effective payload mass of {res.x[0]} [kg]")
print(f"Drag equal to {res.fun[0]} [N]\n")

m_prop = 152
m_fuselage = 260
m_payload = 350
m_structure = 100
m_fuel = m_misc_final - m_payload - m_prop - m_fuselage - m_structure
print(f"Fuel mass left = {m_fuel} [kg]")

cj = 10e-5
R = (cj * D)**-1 * m_fuel * V_cruise / 1000

print(f"Achievable range = {R} [km]")