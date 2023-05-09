import matplotlib.pyplot as plt
import numpy as np
## Constants
rho0 = 1.0 * 10 ** -2  # [kg/m3] Zero elevation density of Mars
T0_low = 145  # [K] Lowest zero elevation temperature
T0_high = 284  # [K] Highest zero elevation temperature
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
    return m, m_lst, j_lst, V, r

for m_misc in list(np.arange(1,3010,1)):
    m_end, m_lst, j_lst, V, r = mass_convergence(m_misc)
    if m_end >= m_goal:
        print(f"\nFinal total mass = {m_end}")
        print(f"Effective payload = {m_misc}\n")
        print(f"Total volume = {V} and radius of sphere = {r}")
        break

plt.plot(j_lst, m_lst)
plt.show()

