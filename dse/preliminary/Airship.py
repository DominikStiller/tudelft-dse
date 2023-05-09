import matplotlib.pyplot as plt
import numpy as np

## Constants
rho0= 1.0 *10 ** -2  # [kg/m3] Zero elevation density of Mars
T0_lst = [145, 284] # [K] List of zero elevation temperatures. Format: T_low, T_high
R_M = 191  # [J/K Kg] Specific Gas constant of Martian Air
R_universal = 8.314  # [J/mol] Universal gas constant
y = 1.2 # [-] Ratio of specific heats
g_M = 3.71  # [m/s2] Gravitational acceleration of Mars
M_H2 = 2.02*10**-3  # [kg/mol] Molar mass of hydrogen
sigma_yield = 1240*10**6  # Yield strength kevlar
p_ratio = 1.005

R_H2 = R_universal/M_H2

i=0
p0 = rho0*T0_lst[i]*R_M
print(p0)
rho_H2 = p_ratio*p0/(R_H2*T0_lst[i])
print(rho_H2)
drho = rho0-rho_H2
m = 3000
m_lst = []
j_lst = []

for j in range(20):
    j_lst.append(j)
    m_lst.append(m)
    V = m/drho
    r = ((3)/(4*np.pi)*V)**(1/3)
    t_min = (p_ratio - 1) * p0 * r / (2 * sigma_yield)
    print(t_min)
    S = 4 * np.pi * r**2
    t = 0.000101
    rho_kevlar = 1380 # [kg/m3]
    m = S*t*rho_kevlar
print(m)
plt.plot(j_lst, m_lst)
plt.show()

