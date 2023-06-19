#!/usr/bin/env python
"""
DSE Group 15

Life Support System (LSS) sizing

Author: Joachim Bron
Created: 25-05-2023
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
# FUNCTIONS
# =============================================================================
def solve(fn, **kwargs):
    """Solve a certain function for the missing argument

    Args:
        fn: function to be solved (needs to be written in the form fn = 0)
        **kwargs: Parameters needed by the function

    Returns:
        The missing parameter value. 

    """
    signature = inspect.signature(fn)
    default_args = {
        k
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }
    arg_to_solve = signature.parameters.keys() - kwargs.keys() - default_args # Find arguments that have been passed
    assert len(arg_to_solve) == 1
    
    def fn_solvable(x):
        args = kwargs | {list(arg_to_solve)[0]: x}
        return fn(**args)
    
    return root(fn_solvable, x0=0.001).x[0]

def print_LSS(lst, units, title):
    """Function used to print a lst in tabulated form 

    Args:
        lst: lst of values and keys
        units: units of the lst
        title: title of the table

    Returns:
        print of attributes of LSS in tabulated form
        
    """
    keys = list(lst.keys())
    values = list(lst.values()) 
    column_titles = ['Variable', 'Value', 'Unit']
    print_lst = [column_titles]
    for i in range(len(keys)):
        local_lst = []
        if not 'object' in str(values[i]):
            local_lst.append(keys[i])
            value = values[i]
            if type(value) != str:
                if value > 1000:
                    value = '%.2E' % Decimal(value)
                else:
                    value = '{0:.3f}'.format(value)
            local_lst.append(value)
            print_lst.append(local_lst)
    for i in range(len(units)):
        print_lst[i+1].append(units[i])

    print(title)
    print(tabulate(print_lst))
    return print_lst
    
def calc_length_tank(V, d):
    """Function used to calculate a parameter of a cylindrical tank with 
    spherical end caps. If used in combination with the "solve" function and l 
    and d are given, V is returned. If V and d are given, l is returned, etc. 

    Args:
        V: tank volume [m3] 
        l: tank length [m]
        d: tank diameter [m]

    Returns:
        The missing tank volume parameter

    """
    return 4 * V / (np.pi * d**2) + 1/3 * d 

def idealgas(p, rho, T, R): 
    """Function used to calculate the missing parameter of the formula of 
    the ideal gas law equation. If used in combination with the "solve" 
    function and p and rho are given, T is returned, etc.
    
    Args:
        p: tank gas pressure [Pa] 
        rho: tank gas density [kg/m3]
        T: tank gas temperature [K]
        R: tank gas constant [J/kg/K]
        
    Returns:
        The missing tank volume parameter

    """
    return p - rho * R * T

def gas_constant(M, R_universal = 8.314): 
    """Calculate the gas constant of a gas
    
    Args:
        M: molar mass of gas [g/mol] 
        
    Returns:
        The gas constant for the specific gas

    """
    return R_universal / M * 1000

def calc_thin_cyl_mass(l, d, t1, t2, rho_material):
    """Function used to calculate the mass of a thin walled cylindrical shape
    
    Args:
        l: tank length [m]
        d: tank diameter [m]
        t1: tank thickness cylindrical part [m]
        t2: tank thickness of spherical end caps [m]
        rho_material: material density [kg/m3]

    Returns:
        The mass of the cylindrical tank

    """
    return rho_material * (np.pi * d * (l-d) * t1 + 4 * np.pi * (d/2)**2 * t2)


def hoop_stress_cyl(sigma_h, p, d, t): 
    """Function used to calculate the missing parameter of the formula of 
    the cylinder hoop stress equation. Use in combination with the 'solve' 
    function.
    
    Args:
        sigma_h: hoop stress [Pa] 
        p: cylinder gauge pressure [Pa]
        d: cylinder diameter [m]
        t: tank gas constant [m]
        
    Returns:
        The missing hoop stress equation parameter when combined with 'solve'

    """
    return sigma_h * t - p * (d / 2)

def longit_stress_cyl(sigma_l, p, d, t): 
    """Function used to calculate the missing parameter of the formula of 
    the cylinder longitudinal stress equation. Use in combination with the 'solve' 
    function.
    
    Args:
        sigma_l: longitudinal stress [Pa] 
        p: cylinder gauge pressure [Pa]
        d: cylinder diameter [m]
        t: tank gas constant [m]
        
    Returns:
        The missing longitudinal stress equation parameter when combined with 'solve'

    """
    return sigma_l * t - p * (d / 2) / 2

def hoop_stress_sph(sigma_h, p, d, t): 
    """Function used to calculate the missing parameter of the formula of 
    the sphere hoop stress equation. Use in combination with the 'solve' 
    function.
    
    Args:
        sigma_h: hoop stress [Pa] 
        p: cylinder gauge pressure [Pa]
        d: cylinder diameter [m]
        t: tank gas constant [m]
        
    Returns:
        The missing sphere hoop stress equation parameter when combined with 'solve'

    """
    return sigma_h * t - p * (d / 2) / 2 

def strain_stress(strain, stress, E):
    """Function used to calculate the missing parameter of the formula of 
    the strain-stress equation. Use in combination with the 'solve' 
    function.
    
    Args:
        strain: strain [-] 
        stress: stress [Pa]
        E: Young's modulus [Pa]
        
    Returns:
        The missing strain-stress equation parameter when combined with 'solve'

    """
    return stress - strain * E


def pressure_liquid(rho, h, g = 3.71): 
    """Function used to calculate the missing parameter of the formula of 
    the sphere hoop stress equation. Use in combination with the 'solve' 
    function.
    
    Args:
        rho: liquid density [kg/m2] 
        h: height of liquid [m]
        g: gravitational acceleration [m/s2]
        
    Returns:
        The missing sphere hoop stress equation parameter when combined with 'solve'

    """
    return rho * g * h

class Astronaut:
    def __init__(self, required_rate_air, required_rate_water, suit_pressure):
        self.required_rate_air = required_rate_air
        self.required_rate_water = required_rate_water
        self.suit_pressure = suit_pressure

class Tank:
    def __init__(self, 
                 content, 
                 content_temperature, 
                 content_pressure, 
                 content_mass,
                 diameter,
                 material,
                 shape = 'cylindrical',
                 SF = 1.1):
        
        # ATTRIBUTES
        # Misc
        self.SF = SF
        
        # Contents
        self.content = content
        self.content_temperature = content_temperature
        self.content_pressure = content_pressure * self.SF
        self.content_mass = content_mass
        
        # Shape and geometric characteristics
        self.diameter = diameter
        self.shape = shape
        
        # Tank
        self.material = material['name']
        self.material_rho = material['rho']
        self.material_sigma_y = material['sigma_y']
        self.material_E = material['E']
        self.material_Kic = material['Kic']
        
    def _calculate_properties_after_setting_subclass(self):
        self.d_over_t = self.diameter / self.t_cyl
        self.thinwalled = str(self.d_over_t > 20)
        self.max_allow_crack = 1 / np.pi * (self.material_Kic / self.material_sigma_y) ** 2
        self.LBB = str(self.max_allow_crack > self.t_cyl)
        self.SF_LBB = self.material_Kic / (self.content_pressure * self.diameter/2 / self.t_cyl / 2 * np.sqrt(np.pi * self.t_cyl))
    
        
class GasTank(Tank):
    def __init__(self, *, gas_M, **kwargs):
        super().__init__(**kwargs)

        # ATTRIBUTES
        self.gas_M = gas_M
        self.gas_R = gas_constant(self.gas_M)
        self.content_rho = solve(idealgas, p = self.content_pressure, T = self.content_temperature, R = self.gas_R)
        self.volume = self.content_mass / self.content_rho 
        self.length = calc_length_tank(V =self.volume, d = self.diameter)
        
        self.t_cyl = solve(hoop_stress_cyl, sigma_h = self.material_sigma_y, 
                                         p = self.content_pressure, 
                                         d = self.diameter)
        self.t_sph = solve(hoop_stress_sph, sigma_h = self.material_sigma_y, 
                                         p = self.content_pressure, 
                                         d = self.diameter)
        
        
        self.empty_mass = calc_thin_cyl_mass( l = self.length, 
                                        d = self.diameter, 
                                        t1 = self.t_cyl, 
                                        t2 = self.t_sph, 
                                        rho_material = self.material_rho)
    

        self._calculate_properties_after_setting_subclass()
        
class LiquidTank(Tank):
    def __init__(self, *, content_rho, **kwargs):
        super().__init__(**kwargs)

        # ATTRIBUTES
        self.content_rho = content_rho
        self.volume = self.content_mass / self.content_rho 
        self.length = calc_length_tank(V =self.volume, d = self.diameter)
        
        self.content_pressure = pressure_liquid(self.content_rho, self.length) * self.SF
        self.t_cyl = max(solve(hoop_stress_cyl, sigma_h = self.material_sigma_y, 
                                         p = self.content_pressure, 
                                         d = self.diameter), 2e-3)
        self.t_sph = self.t_cyl # IS THIS VALID?
        
        self.empty_mass = calc_thin_cyl_mass( l = self.length, 
                                        d = self.diameter, 
                                        t1 = self.t_cyl, 
                                        t2 = self.t_sph, 
                                        rho_material = self.material_rho)
        
        self._calculate_properties_after_setting_subclass()

        
class LifeSupportSystem:
    def __init__(self, air_leakage_rate, air_days_needed, water_days_needed, num_tanks_air, num_tanks_water, num_astronauts):
        
        self.air_leakage_rate = air_leakage_rate
        self.air_days_needed = air_days_needed
        self.water_days_needed = water_days_needed
        self.num_tanks_air = num_tanks_air
        self.num_tanks_water = num_tanks_water
        self.num_astronauts = num_astronauts
        
        self.Astronauts = Astronaut(required_rate_air, required_rate_water, suit_pressure)
        self.required_rate_air_LSS = (self.Astronauts.required_rate_air + self.air_leakage_rate) * self.num_astronauts
        self.required_rate_water_LSS = self.Astronauts.required_rate_water * self.num_astronauts
        
        self.required_mass_air = self.required_rate_air_LSS * self.air_days_needed
        self.required_mass_water = self.required_rate_water_LSS * self.water_days_needed
        
        self.mass_per_tank_O2 = self.required_mass_air / self.num_tanks_air
        self.mass_per_tank_H20 = self.required_mass_water / self.num_tanks_water
        
        self.O2Tank = GasTank(content = 'O2', 
                              content_temperature = tank_temperature_air, 
                              content_pressure = tank_pressure, 
                              content_mass = self.mass_per_tank_O2, 
                              diameter = tank_diameter,
                              material = Steel,
                              SF = 1,
                              gas_M = O2_M)
        
        self.H2OTank = LiquidTank(content = 'H2O', 
                                  content_temperature = tank_temperature_water, 
                                  content_pressure = tank_pressure, 
                                  content_mass = self.mass_per_tank_H20, 
                                  diameter = tank_diameter,
                                  material = Steel,
                                  SF = 1.1,
                                  content_rho = density_water)
        self.total_mass = self.required_mass_air + self.required_mass_water + self.O2Tank.empty_mass * self.num_tanks_air + self.H2OTank.empty_mass * self.num_tanks_water
        self.empty_mass = self.O2Tank.empty_mass * self.num_tanks_air + self.H2OTank.empty_mass * self.num_tanks_water

# =============================================================================
# INPUT PARAMETERS
# =============================================================================
# Fluids properties
O2_M =                  32      #[g/mol]
Mars_atm_pressure =     6.518e2 #[Pa] #Unusued atm. Should be used to get gauge pressure by subtracting from internal tank pressure, but is negligibly small
density_water =         997     #[kg/m3]

# Astronauts
air_leakage_rate =      40e-3   #[kg/day/module]
required_rate_air =     840e-3  #[kg/day/module]
required_rate_water =   4e-3 * density_water #[kg/day]
air_days_needed =       1.5       #[day]
water_days_needed =     1       #[day]
num_astronauts =        2       #[-]
suit_pressure =         20e3    #[Pa]

# Tanks
num_tanks_air =             1           #[-]
num_tanks_water =           1           #[-]
tank_temperature_air =      18+273.15   #[K]
tank_temperature_water =    18+273.15   #[K]
tank_pressure =             20000e3     #[Pa]
tank_diameter =             14e-2       #[m]

# Material properties
steel_rho = 8000        #[kg/m3]
steel_E = 193e9         #[Pa]
steel_sigma_y = 215e6   #[Pa]
steel_Kic = 100e6       #[Pa m^1/2]

# =============================================================================
# Review these values 
# =============================================================================
cfrp_rho = 2000         #[kg/m3]
cfrp_E = 220e9            #[Pa]
cfrp_sigma_y = 0.3 * 880e6    #[Pa]
cfrp_Kic = 10e6   #[Pa m^1/2]
# =============================================================================

Steel =  {'name': 'Steel', 'rho': steel_rho, 'sigma_y': steel_sigma_y, 'E': steel_E, 'Kic': steel_Kic}
CFRP =  {'name': 'CFRP', 'rho': cfrp_rho, 'sigma_y': cfrp_sigma_y, 'E': cfrp_E, 'Kic': cfrp_Kic}

if __name__ == '__main__':
    

        
    LSS = LifeSupportSystem(air_leakage_rate, air_days_needed, water_days_needed, num_tanks_air, num_tanks_water, num_astronauts)
    # pprint(inspect.getmembers(LSS))
    # pprint(inspect.getmembers(LSS.O2Tank))
    # pprint(inspect.getmembers(LSS.H20Tank))
    # print(tabulate(inspect.getmembers(LSS)[3]))
    
    def convert_to_lst(lst_in):
        lst = inspect.getmembers(lst_in, lambda a:not(inspect.isroutine(a)))
        lst = dict([a for a in lst if not(a[0].startswith('__') and a[0].endswith('__'))])
        return lst
    
    lst1 = convert_to_lst(LSS)
    lst2 = convert_to_lst(LSS.Astronauts)
    lst3 = convert_to_lst(LSS.O2Tank)
    lst4 = convert_to_lst(LSS.H2OTank)
    
    units1 = ['-', 'kg/day', 'kg', 'kg', 'kg', '-', '-', '-', 'kg', 'kg', 'kg/day', 'l/day', 'kg', '-']
    units2 = ['kg/day', 'l/day', 'Pa']
    units3 = ['-', '-', '-', '-', 'kg', 'Pa', 'kg/m3', 'K', '-', 'm', 'kg', 'g/mol', 'J/kg/K', 'm', '-', 'Pa', 'Pa m1/2', 'kg/m3', 'Pa', 'm', '-', 'm', 'm', '-', 'm3']
    units4 = ['-', '-', '-', '-', 'kg', 'Pa', 'kg/m3', 'K', '-', 'm', 'kg', 'm', '-', 'Pa', 'Pa m1/2', 'kg/m3', 'Pa', 'm', '-', 'm', 'm', '-', 'm3']

    print_LSS(lst1, units1, 'LSS')
    print_LSS(lst2, units2, 'Astronauts')
    table = print_LSS(lst3, units3, 'O2 tank')
    print_LSS(lst4, units4, 'H2O tank')
    
    
    import matplotlib.pyplot as plt
    
    l_over_d = 3
    rho_wall = steel_rho
    R = 260
    T = 290
    n_tanks = 3
    C1 = 1 / R /T
    
    C2 = 1/4 * np.pi * (l_over_d - 1/3)
    
    C3 = 1 / 4 / steel_sigma_y
    m_o2 = 1.76 / n_tanks
    
    d_min = 5e-2
    d_max = 30e-2
    d_range = np.linspace(d_min,d_max,100)
    d_range = np.reshape(d_range, (-1,1))
    
    t_min = 1e-3
    t_max = 15e-3
    t_range = np.linspace(t_min,t_max,100)
    
    m_wall = np.pi * d_range ** 2 * t_range * (2 * l_over_d - 1) * rho_wall
    
    xx, yy = np.meshgrid(t_range,d_range)
    c1 = 20 * t_range
    c2 = (m_o2 * C3 / (C2 * C1 * t_range)) ** (0.5)
    c5 = (m_o2 * R * T / (steel_sigma_y * np.pi * (l_over_d - 1/3) * t_range)) ** 0.5
    #c4 constraint is due to crack LBB
    # c4 = 1 / np.pi * (steel_Kic / steel_sigma_y) ** 2
    # plt.imshow(m_wall,extent=(t_min*1000,t_max*1000,100*d_min,100*d_max),aspect='auto')
    plt.contourf(xx*1000,yy*1000,m_wall, levels = 20)

    plt.colorbar()
    plt.plot(t_range*1000,c1*1000)
    plt.plot(t_range*1000,c1*1000)
    plt.fill_between(t_range*1000, 0, np.maximum(c1*1000, c5*1000), alpha=1, color="white")
    plt.plot(t_range*1000,c5*1000)
    cs = plt.contour(xx*1000,yy*1000,m_wall, levels = 20)
    plt.clabel(cs)
    # plt.axvline(c4 * 1000)
    # c3_t_range = np.repeat(t_min, 100)
    # c3 = d_range
    
    m_wall_min = np.pi * c2 ** 2 * t_range * (2 * l_over_d - 1) * rho_wall
    print(m_wall_min)
    print(m_wall_min * n_tanks)
    
