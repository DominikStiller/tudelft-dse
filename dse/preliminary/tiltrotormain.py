from RotorEngineSizing import RadiusMassElementMomentum
from AircraftEstimating import Class2Weight
from cruise_sizing import max_tipSpeed

#Top Level Requirement:
Mass_design = 3000 #kg
V_cr = 112 #m/s
Mass_payload = 400 #kg
#Sizing of Rotor:
margin = 0.1
Mass_thrust = Mass_design*(1+margin) #kg
TotalRotors = 4
coaxial = True
BladePerRotor = 6

Mass_solar = 0 #kg

#Weight Estimation Inputs
N_ult = 4.4
AR = 20
wingbraced = True

#Constants for Battery Sizing
E_density = 333.33 #W/kg
P_density_TO = 700 #Wh/kg Power Density of the devices used for take-off
E_density_TO =  300 #W/kg
Volume_bat = 450 # Wh/L

#Set up Calculations
TipSpeed = -0.88*V_cr + 268.87

R, T, Horsepower, Power, mass_rotors = RadiusMassElementMomentum(Mass_thrust, TotalRotors, BladePerRotor, coaxial, TipSpeed)
Class2Weight(R, mass_rotors, Mass_design, N_ult, AR, wingbraced, V_cr, E_density, P_density_TO, E_density_TO, Mass_payload, Mass_solar)

