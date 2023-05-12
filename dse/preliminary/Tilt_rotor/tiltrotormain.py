from PerformanceAnalysis import PayloadRange, AircraftClimbPerf, RotorClimbPerf
from AircraftEstimating import Class2Weight, DragEstimation
from RotorEngineSizing import RadiusMassElementMomentum
from power_sizing import size_power_subsystem
from cruise_sizing import area_and_thrust
from constants import const


if __name__ == '__main__':
    # Sizing of Rotor:
    margin = 0.1

    Mass_design = 3000  # kg
    Mass_thrust = Mass_design*(1+margin)  # kg
    Mass_payload = 400  # kg

    DesignRange = 1000 #km
    TakeoffTime = 10/60 * 1

    TotalRotors = 4
    coaxial = True
    BladePerRotor = 6

    N_ult = 4.4
    AR = 20
    wingbraced = True
    V_cr = 400/3.6  # m/s
    E_density = 333.33  # W/kg
    P_density_cr = 300  # wh/kg
    P_density_TO = 700  # Wh/kg Power Density of the devices used for take-off
    E_density_TO = 300  # W/kg
    Volume_bat = 450  # Wh/L
    # Set up Calculations
    TipSpeed = -0.88*V_cr + 268.87
    diff=100
    while diff>0.01:
        Mass_thrust = Mass_design * (1 + margin)  # kg
        R, takeOffThrustPerEngine, Horsepower, Power, mass_rotors = RadiusMassElementMomentum(Mass_thrust, TotalRotors, BladePerRotor, coaxial, TipSpeed, print_results=True)

        # Calculate wing area
        S, cruiseThrust = area_and_thrust(0, const['cl'], const['cd'], Mass_design, 0.5*const['airDensity']*V_cr**2)

        # Size batteries for cruise and solar panels
        cruiseThrust += DragEstimation(R=R, Swing=S, t2c=const['t/c'], Vcr=V_cr, visc_cr=const['visc_cr'], AR=AR)
        batteryMass, panelMass, powerSurplus = size_power_subsystem(R, takeOffThrustPerEngine, cruiseThrust, DesignRange / V_cr + TakeoffTime, TakeoffTime, S,
                                                                      plot=True)

        # Calculate weights
        Range, wingWeight, tailWeight, bodyWeight = Class2Weight(R, mass_rotors, Mass_design, N_ult, AR,
                                                                                        wingbraced, V_cr, E_density, P_density_TO,
                                                                                        E_density_TO, Mass_payload, panelMass)

        # Payload-Range
        PayloadRange(R, mass_rotors, Mass_design, N_ult, AR, wingbraced, V_cr, E_density, P_density_TO, E_density_TO, panelMass, Mass_payload, DesignRange)

        # Climb Performance
        ROC_cruise = AircraftClimbPerf(batteryMass, P_density_cr, Mass_design, R, V_cr)
        ROC_rotor = RotorClimbPerf(Mass_design, R, TotalRotors)
        Totalweight = mass_rotors+panelMass+batteryMass+wingWeight+tailWeight+bodyWeight
        diff = abs((Totalweight-Mass_design)/Mass_design)
        print(Totalweight)
        Mass_design = Totalweight
        #diff = 0.001
