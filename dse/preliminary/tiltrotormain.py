from RotorEngineSizing import RadiusMassElementMomentum
from AircraftEstimating import Class2Weight
from cruise_sizing import area_and_thrust
from constants import const
from dse.preliminary.PerformanceAnalysis import PayloadRange

if __name__ == '__main__':
    # Sizing of Rotor:
    margin = 0.1
    Mass_thrust = 3000*(1+margin)  # kg
    Mass_design = 3000  # kg
    Mass_payload = 400  # kg
    Mass_solar = 0#kg

    TotalRotors = 4
    coaxial = True
    TipSpeed = 160  # m/s
    BladePerRotor = 6

    N_ult = 4.4
    AR = 20
    wingbraced = True
    V_cr = 112  # m/s
    E_density = 333.33  # W/kg
    P_density_TO = 700  # Wh/kg Power Density of the devices used for take-off
    E_density_TO = 300  # W/kg
    Volume_bat = 450  # Wh/L

#Set up Calculations
    TipSpeed = -0.88*V_cr + 268.87
    R, T, Horsepower, Power, mass_rotors = RadiusMassElementMomentum(Mass_thrust, TotalRotors, BladePerRotor, coaxial, TipSpeed, print_results=True)
    Class2Weight(R, mass_rotors, Mass_design, N_ult, AR, wingbraced, V_cr, E_density, P_density_TO, E_density_TO, Mass_payload, Mass_solar)

    # Calculate wing area
    S, cruiseThrust = area_and_thrust(0, const['cl'], const['cd'], Mass_design, 0.5*const['airDensity']*V_cr**2)
    Mass_solar = const['solarPanelDensity'] * S

    # Calculate weights
    Range, wingWeight, tailWeight, bodyWeight, controlSurfacesWeight = Class2Weight(R, mass_rotors, Mass_design, N_ult, AR,
                                                                             wingbraced, V_cr, E_density, P_density_TO,
                                                                             E_density_TO, Mass_payload, Mass_solar)

    #Payload-Range
    PayloadRange(R, mass_rotors, Mass_design, N_ult, AR, wingbraced, V_cr, E_density, P_density_TO, E_density_TO,
                 Mass_solar, Mass_payload, minRange=1000)

