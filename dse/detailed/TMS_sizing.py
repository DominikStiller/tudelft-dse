#!/usr/bin/env python
"""
DSE Group 15

Thermal Management System (TMS) sizing

Author: Joachim Bron
Created: 12-06-2023
"""
# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import inspect
from functools import partial
from scipy.optimize import root
from pprint import pprint
import inspect
from tabulate import tabulate
from decimal import Decimal

# =============================================================================
# Functions
# =============================================================================
def calc_heat_in_solar(Is, alpha, Fsr, Aperp):
    return Is * alpha * Fsr * Aperp

def calc_heat_in_albedo(b, Fsr, Aperp, Qsun):
    return b * Qsun

def calc_heat_in_infrared(Aperp, Jp):
    return Jp * Aperp

def calc_heat_int_engine(Pout, eta):
    return Pout * (1-eta)/eta

def calc_heat_int_battery(Pout, eta):
    return Pout * (1-eta)/eta

def calc_heat_out_convection(h, Tac, Tamb, Aconv):
    return h * (Tac - Tamb) * Aconv

def calc_heat_out_radiation(epsilon_ac, Tac, Tamb, Arad, sigma = 5.6704e-8):
    return sigma * epsilon_ac * (Tac ** 4 - Tamb ** 4) * Arad

def calc_eqm_temperature(Q_in_solar, Q_in_albedo, Q_in_infrared, Q_in_int, h, Aconv, Tamb, epsilon_ac, Arad, sigma = 5.6704e-8):
    a = sigma * epsilon_ac * Arad
    b = 0
    c = 0 
    d = h * Aconv
    e = - (Q_in_solar + Q_in_albedo + Q_in_infrared + Q_in_int + h * Aconv * Tamb + sigma * epsilon_ac * Arad * Tamb ** 4)

    p = [a, b, c, d, e]
    roots = np.roots(p)
    return roots

def calc_Cp_func_T(T):
    if T > 1000 or T < 200:
        print("Warning: Cp function not valid for T outside of [200,1000]")
    a1 = 4.452e2
    a2 = 1.697
    a3 = -1.346e-3
    a4 = 4.646e-7
    a5 = -2.715e-11
    return a1 + a2 * T + a3 * T ** 2 + a4 * T ** 3 + a5 * T **4

def calc_kf_func_T(T):
    if T > 1670 or T < 185:
        print("Warning: kf function not valid for T outside of [185,1670]")
    a1 = -7.215e-3
    a2 = 8.015e-5
    a3 = 5.477e-9
    a4 = -1.053e-11
    return a1 + a2 * T + a3 * T ** 2 + a4 * T ** 3 

def calc_Prandtl_number(Cp, kf, mu):
    return Cp * mu / kf

def calc_Reynolds_number(rho, v, L, mu):
    return rho * v * L / mu

def calc_Nusselt_number(Re, Pr):
    return 0.3 + (0.62 * Re ** (1/2) * Pr ** (1/3))/(1+(0.4/Pr)**(2/3))**(1/4) * (1+(Re/282000)**(5/8))**(4/5)

def calc_h(Nu, L, kf):
    return Nu * kf / L

#%%
# =============================================================================
# Batteries, worst case [hot, cold] 
# =============================================================================

# Inputs 

# Boundary conditions
Tamb = [270.730, 183.231]  #[K]
Aw = 134 #[m2]
Ah = np.nan #[m2] # Not of interest
Av = Ah/2 #[m2]
Af = (4 * 1.5) * 3 + (1.5 * 1.5) * 2 #[m2] #l_f * w_f * 3 sides + w_f * h_f * 2 sides

# Surface area battery: length of 2 m, quarter circle of radius 0.5 m, 2m x (2 pi 0.5 / 2 + 4 * 0.5 ) + pi 0.5^2 = 
V_bat = 0.5 #[m3]
l_bat = 2 #[m]
A_bat = V_bat/l_bat #[m2]
r_bat = np.sqrt(2 * A_bat / np.pi) #[m]

Abat_inside = 2 * (2 * l_bat * r_bat) #[m2] 2 sides per battery, 2 batteries
Abat_outside = l_bat * (2 * np.pi * r_bat) / 2 #[m2]
Abat_front_aft = np.pi * r_bat ** 2 
Abat_surface = Abat_inside + Abat_outside + Abat_front_aft 

# Atailboom = 10 * np.pi * 0.2 # unused

# Engine surface area
r_eng = 0.1 #[m]
l_eng = 1.99 #[m] same as tip chord
Aeng = 2 * (2 * np.pi * r_eng * l_eng + 2 * np.pi * r_eng ** 2)

fw = 0
fbi = 14 # number of areas inside if battery is cut up
fbo = 5 * 2 + 1 # *2 because each radiator plate has two sides 
fbaf = 1

rho_material_radiator = 2700
t_radiator_cells = 1e-3
m_radiator_battery = (fbo - 1) / 2 * Abat_outside * t_radiator_cells * rho_material_radiator


Arad_out = [fbo * Abat_outside + fbi * Abat_inside + fw * Aw, fbo * Abat_outside + fbi * Abat_inside + fw * Aw]
Aconv = [fbi * Abat_inside + fbo * Abat_outside + fbaf * Abat_front_aft + fw * Aw, 0] #fbi * Abat_inside + fbo * Abat_outside + fbaf * Abat_front_aft + fw * Aw] 
Arad_in = 0 # + fw * Aw #+ Aw # Same as A_perp in paper

# Air parameters
rho = 1.482e-2
mu = 9.82e-6
v = [2,111] #m/s
L = 2

Cp = [calc_Cp_func_T(Tamb[0]), calc_Cp_func_T(Tamb[1])]
kf = [calc_kf_func_T(Tamb[0]), calc_kf_func_T(Tamb[1]) ] 
Pr = [calc_Prandtl_number(Cp[0],kf[0],mu), calc_Prandtl_number(Cp[1],kf[1],mu)]
Re = [calc_Reynolds_number(rho, v[0], L, mu), calc_Reynolds_number(rho, v[1], L, mu)]
Nu = [calc_Nusselt_number(Re[0], Pr[0]), calc_Nusselt_number(Re[1], Pr[1])]
h = [calc_h(Nu[0], L, kf[0]), calc_h(Nu[1], L, kf[1])]

# Direct solar radiation
Is = [587.424, 0] #[W/m2]
alpha = 0.2 #[-] #change 
Fsr = 1 #[-] #assumed

# Indirect solar radiation
b = 0.4 #[-] #const
Fsr = 1 #[-] #assumed

# Planetary radiation
Jp = [70.72, 16.881] #[W/m2] #const

# No conduction was sized

# Electronics
Q_in_int_electronics = 25 + 25 #[W] # 25 -> flight computer + 25 assumed for the rest of electrical

# Engine
Pout_eng_to = 312e3 #[W] # Power output engine take-off
Pout_eng_cruise = 265e3 #[W] # Power output engine cruise
eta_eng_to = 0.98 ** 3 * 0.9 #[-] # Gearbox and electric engine
eta_eng_cruise = 0.98 ** 3 * 0.93 #[-] # Gearbox and electric engine

# Batteries 
Pout_bat_to = Pout_eng_to / eta_eng_to #[W] Pout_bat = Pin_eng
Pout_bat_cruise = Pout_eng_cruise / eta_eng_cruise #[W]
eta_bat = 0.95 #[-]

# Heat output
# Convection radiation
h = h #[W/m2/K]
Aconv = Aconv

# Thermal radiation
Fsr = Fsr #[-]
epsilon_ac = 0.85 #[-]
Arad_out = Arad_out

# Heat from engines and batteries (applies to both hot and cold scenario)
# Q_in_int_engine_to = calc_heat_int_engine(Pout_eng_to, eta_eng_to) #[W]
# Q_in_int_engine_cruise = calc_heat_int_engine(Pout_eng_cruise, eta_eng_cruise) #[W]
Q_in_int_battery_to = calc_heat_int_battery(Pout_bat_to, eta_bat) #[W]
Q_in_int_battery_cruise = calc_heat_int_battery(Pout_bat_cruise, eta_bat) #[W]
Q_in_int_electronics = Q_in_int_electronics #[W]


# Worst case hot
Q_in_solar_h = calc_heat_in_solar(Is[0], alpha, Fsr, Arad_in) #[W]
Q_in_albedo_h = calc_heat_in_albedo(b, Fsr, Arad_in, Q_in_solar_h) #[W]
Q_in_infrared_h = calc_heat_in_infrared(Arad_in, Jp[0]) #[W]
Q_in_int_total_h = max(Q_in_int_battery_to, Q_in_int_battery_cruise) + Q_in_int_electronics  #[W]

Tac_h = calc_eqm_temperature(Q_in_solar = Q_in_solar_h, 
                           Q_in_albedo = Q_in_albedo_h, 
                           Q_in_infrared = Q_in_infrared_h, 
                           Q_in_int = Q_in_int_total_h, 
                           h = h[0], Aconv = Aconv[0], Tamb = Tamb[0], 
                           epsilon_ac = epsilon_ac, Arad = Arad_out[0]) #[K], [deg C]

Q_out_conv_h = calc_heat_out_convection(h = h[0], Tac = Tac_h[-1].real, Tamb = Tamb[0], Aconv = Aconv[0])
Q_out_rad_h = calc_heat_out_radiation(epsilon_ac = epsilon_ac, Tac = Tac_h[-1].real, Tamb = Tamb[0], Arad = Arad_out[0])


# Worst case cold
Q_in_solar_c = calc_heat_in_solar(Is[1], alpha, Fsr, Arad_in) #[W]
Q_in_albedo_c = calc_heat_in_albedo(b, Fsr, Arad_in, Q_in_solar_c) #[W]
Q_in_infrared_c = calc_heat_in_infrared(Arad_in, Jp[1]) #[W]
Q_in_int_total_c = min(Q_in_int_battery_to, Q_in_int_battery_cruise) + Q_in_int_electronics  #[W]

Tac_c = calc_eqm_temperature(Q_in_solar = Q_in_solar_c, 
                           Q_in_albedo = Q_in_albedo_c, 
                           Q_in_infrared = Q_in_infrared_c, 
                           Q_in_int = Q_in_int_total_c, 
                           h = h[1], Aconv = Aconv[1], Tamb = Tamb[1], 
                           epsilon_ac = epsilon_ac, Arad = Arad_out[1]) #[K], [deg C]

Q_out_conv_c = calc_heat_out_convection(h = h[1], Tac = Tac_h[-1].real, Tamb = Tamb[1], Aconv = Aconv[1])
Q_out_rad_c = calc_heat_out_radiation(epsilon_ac = epsilon_ac, Tac = Tac_h[-1].real, Tamb = Tamb[1], Arad = Arad_out[1])


print(f"\n {Tac_h - 273} \n")
print(f"\n {Tac_c - 273} \n")


#%%
# =============================================================================
# Engines, worst case [hot, cold] 
# =============================================================================

# Inputs

# Boundary conditions
Tamb = [270.730, 183.231]  #[K]
Aw = 134 #[m2]
Ah = np.nan #[m2] # Not of interest
Av = Ah/2 #[m2]
Af = (4 * 1.5) * 3 + (1.5 * 1.5) * 2 #[m2] #l_f * w_f * 3 sides + w_f * h_f * 2 sides

# Surface area battery: length of 2 m, quarter circle of radius 0.5 m, 2m x (2 pi 0.5 / 2 + 4 * 0.5 ) + pi 0.5^2 = 
V_bat = 0.5 #[m3]
l_bat = 2 #[m]
A_bat = V_bat/l_bat #[m2]
r_bat = np.sqrt(2 * A_bat / np.pi) #[m]

Abat_inside = 2 * (2 * l_bat * r_bat) #[m2] 2 sides per battery, 2 batteries
Abat_outside = l_bat * (2 * np.pi * r_bat) / 2 #[m2]
Abat_front_aft = np.pi * r_bat ** 2 
Abat_surface = Abat_inside + Abat_outside + Abat_front_aft 

# Atailboom = 10 * np.pi * 0.2 # unused

# Engine surface area
r_eng = 0.1 #[m]
l_eng = 1.99 #[m] same as tip chord
Aeng = 2 * (2 * np.pi * r_eng * l_eng + 2 * np.pi * r_eng ** 2)

# Engine radiator
r_radiator = 0.2
l_radiator = l_eng
n_radial = 30 # number of radial plates
w_avg_cell = 2 * np.pi / n_radial * ((r_radiator - r_eng)/2 + r_eng)

A_radial = n_radial * (r_radiator - r_eng) * l_radiator

n_circ = 10
A_circ = 0
for x in range(n_circ):
    A_circ += 2 * np.pi * (x * (r_radiator - r_eng) / n_circ + r_eng) 

h_avg_cell = (r_radiator - r_eng) / n_circ
Aradiator = 2 * (A_radial + A_circ) * 2 #The last *2 accounts for both sides of a cell wall cooling


rho_material_radiator = 2700 #aluminium
t_radiator_cells = 5e-4
m_radiator = Aradiator * t_radiator_cells * rho_material_radiator

fe = 1 
fw = 0
Arad_out = Aeng + Aradiator + fw * Aw
Aconv = [fe * Aeng + Aradiator +fw * Aw, fe * Aeng + Aradiator + fw * Aw] 
Arad_in = Aeng #+ fw * Aw

# Air parameters
rho = 1.482e-2
mu = 9.82e-6
v = [2,111] #m/s
L = 2

Cp = [calc_Cp_func_T(Tamb[0]), calc_Cp_func_T(Tamb[1])]
kf = [calc_kf_func_T(Tamb[0]), calc_kf_func_T(Tamb[1]) ] 
Pr = [calc_Prandtl_number(Cp[0],kf[0],mu), calc_Prandtl_number(Cp[1],kf[1],mu)]
Re = [calc_Reynolds_number(rho, v[0], L, mu), calc_Reynolds_number(rho, v[1], L, mu)]
Nu = [calc_Nusselt_number(Re[0], Pr[0]), calc_Nusselt_number(Re[1], Pr[1])]
h = [calc_h(Nu[0], L, kf[0]), calc_h(Nu[1], L, kf[1])]

# Direct solar radiation
Is = [587.424, 0] #[W/m2]
alpha = 0.2 #[-] #change 
Fsr = 1 #[-] #assumed

# Indirect solar radiation
b = 0.4 #[-] #const
Fsr = 1 #[-] #assumed

# Planetary radiation
Jp = [70.72, 16.881] #[W/m2] #const

# No conduction was sized

# Electronics
Q_in_int_electronics = 25 + 25 #[W] # 25 -> flight computer + 25 assumed for the rest of electrical

# Engine
Pout_eng_to = 312e3 #[W] # Power output engine take-off
Pout_eng_cruise = 265e3 #[W] # Power output engine cruise
eta_eng_to = 0.98 ** 3 * 0.9 #[-] # Gearbox and electric engine
eta_eng_cruise = 0.98 ** 3 * 0.93 #[-] # Gearbox and electric engine

# Batteries 
Pout_bat_to = Pout_eng_to / eta_eng_to #[W] Pout_bat = Pin_eng
Pout_bat_cruise = Pout_eng_cruise / eta_eng_cruise #[W]
eta_bat = 0.95 #[-]

# Heat output
# Convection radiation
h = h #[W/m2/K]
Aconv = Aconv

# Thermal radiation
Fsr = Fsr #[-]
epsilon_ac = 0.85 #[-]
Arad_out = Arad_out


# Heat from engines and batteries (applies to both hot and cold scenario)
Q_in_int_engine_to = calc_heat_int_engine(Pout_eng_to, eta_eng_to) #[W]
Q_in_int_engine_cruise = calc_heat_int_engine(Pout_eng_cruise, eta_eng_cruise) #[W]
# Q_in_int_battery_to = calc_heat_int_battery(Pout_bat_to, eta_bat) #[W]
# Q_in_int_battery_cruise = calc_heat_int_battery(Pout_bat_cruise, eta_bat) #[W]
# Q_in_int_electronics = Q_in_int_electronics #[W]


# Worst case hot
Q_in_solar_h = calc_heat_in_solar(Is[0], alpha, Fsr, Arad_in) #[W]
Q_in_albedo_h = calc_heat_in_albedo(b, Fsr, Arad_in, Q_in_solar_h) #[W]
Q_in_infrared_h = calc_heat_in_infrared(Arad_in, Jp[0]) #[W]
Q_in_int_total_h = max(Q_in_int_engine_to, Q_in_int_engine_cruise)  #[W]

Tac_h = calc_eqm_temperature(Q_in_solar = Q_in_solar_h, 
                           Q_in_albedo = Q_in_albedo_h, 
                           Q_in_infrared = Q_in_infrared_h, 
                           Q_in_int = Q_in_int_total_h, 
                           h = h[0], Aconv = Aconv[0], Tamb = Tamb[0], 
                           epsilon_ac = epsilon_ac, Arad = Arad_out) #[K], [deg C]

Q_out_conv_h = calc_heat_out_convection(h = h[0], Tac = Tac_h[-1].real, Tamb = Tamb[0], Aconv = Aconv[0])
Q_out_rad_h = calc_heat_out_radiation(epsilon_ac = epsilon_ac, Tac = Tac_h[-1].real, Tamb = Tamb[0], Arad = Arad_out)


# Worst case cold
Q_in_solar_c = calc_heat_in_solar(Is[1], alpha, Fsr, Arad_in) #[W]
Q_in_albedo_c = calc_heat_in_albedo(b, Fsr, Arad_in, Q_in_solar_c) #[W]
Q_in_infrared_c = calc_heat_in_infrared(Arad_in, Jp[1]) #[W]
Q_in_int_total_c = min(Q_in_int_engine_to, Q_in_int_engine_cruise) #[W]

Tac_c = calc_eqm_temperature(Q_in_solar = Q_in_solar_c, 
                           Q_in_albedo = Q_in_albedo_c, 
                           Q_in_infrared = Q_in_infrared_c, 
                           Q_in_int = Q_in_int_total_c, 
                           h = h[1], Aconv = Aconv[1], Tamb = Tamb[1], 
                           epsilon_ac = epsilon_ac, Arad = Arad_out) #[K], [deg C]

Q_out_conv_c = calc_heat_out_convection(h = h[1], Tac = Tac_h[-1].real, Tamb = Tamb[1], Aconv = Aconv[1])
Q_out_rad_c = calc_heat_out_radiation(epsilon_ac = epsilon_ac, Tac = Tac_h[-1].real, Tamb = Tamb[1], Arad = Arad_out)


print(f"\n {Tac_h - 273} \n")
print(f"\n {Tac_c - 273} \n")

# Pipes mass
span = 42 
pipe_diam = 3e-2
pipe_thickness = 1.5e-3
pipe_crosssectional_area = np.pi * pipe_diam * pipe_thickness
pipe_volume = 1.5 * span * 2 * 2 * pipe_crosssectional_area #1.5 safety factor, 2 * 2 because 2 times half a wing for full wing and 2 times to go to engine and come back to fuselage
pipe_mass = pipe_volume * 2700 #2700 density of alum


if __name__ == '__main__':
    pass