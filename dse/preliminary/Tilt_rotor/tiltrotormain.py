from PerformanceAnalysis import PayloadRange, AircraftClimbPerf, RotorClimbPerf
from AircraftEstimating import Class2Weight, DragEstimation
from RotorEngineSizing import RadiusMassElementMomentum
from constants import const, aircraftParameters
from power_sizing import size_power_subsystem
from cruise_sizing import area_and_thrust
import xarray as xr
import numpy as np


def design(iterate=False):
    diff = 100
    while diff > 0.01:
        designMass = aircraftParameters['totalMass']+0.00000001
        tipSpeed = -0.88 * const['cruiseSpeed'] + 268.87
        Mass_thrust = Mass_design * (1 + const['margin'])  # kg
        aircraftParameters['rotorRadius'], takeOffThrustPerEngine, aircraftParameters['horsepowerPerEngine'], \
        aircraftParameters['totalPower'], aircraftParameters['rotorMass'] = \
            RadiusMassElementMomentum(M=Mass_thrust, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=tipSpeed, print_results=True)

        # Calculate wing area
        aircraftParameters['wingArea'], aircraftParameters['cruiseThrust'] = \
            area_and_thrust(0, const['cl'], const['cd'], Mass_design,
                            0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2)

        aircraftParameters['wingspan'] = 3 * aircraftParameters['rotorRadius']
        aircraftParameters['chord'] = aircraftParameters['wingArea'] / aircraftParameters['wingspan']
        aircraftParameters['AR'] = aircraftParameters['wingspan'] / aircraftParameters['chord']

        if aircraftParameters['wingArea'] / aircraftParameters['wingspan'] < aircraftParameters['ARmin']:
            aircraftParameters['wingspan'] = np.sqrt(aircraftParameters['wingArea'] * aircraftParameters['ARmin'])
            aircraftParameters['chord'] = aircraftParameters['wingArea'] / aircraftParameters['wingspan']

        print(f'Wing area = {aircraftParameters["wingArea"]}[m2]')
        print(f'Wingspan = {aircraftParameters["wingspan"]}[m]')
        print(f'Chord = {aircraftParameters["chord"]}[m]')
        print(f'AR = {aircraftParameters["wingspan"] / aircraftParameters["chord"]}')

        # Size batteries for cruise and solar panels
        aircraftParameters['cruiseThrust'] += DragEstimation(R=aircraftParameters['rotorRadius'],
                                                         Swing=aircraftParameters['wingArea'], t2c=const['t/c'],
                                                         Vcr=const['cruiseSpeed'], visc_cr=const['visc_cr'],
                                                         AR=aircraftParameters['AR'])

        aircraftParameters['batteryMass'], aircraftParameters['panelMass'], powerSurplus = \
            size_power_subsystem(aircraftParameters['rotorRadius'], takeOffThrustPerEngine,
                                 aircraftParameters['cruiseThrust'],
                                 const['designRange'] / const['cruiseSpeed'] + const['takeOffTime'],
                                 const['takeOffTime'], aircraftParameters['wingArea'], plot=True)

        # Calculate weights
        Range, aircraftParameters['wingMass'], aircraftParameters['tailMass'], aircraftParameters['bodyMass'] = \
            Class2Weight(aircraftParameters['rotorRadius'], aircraftParameters['rotorMass'], const['maxMass'],
                         const['ultimateLoad'], aircraftParameters['AR'], aircraftParameters['wingbraced'],
                         const['cruiseSpeed'], const['batteryEnergyDensity'], const['batteryPowerDensity'],
                         const['batteryEnergyDensity'], const['payloadMass'], aircraftParameters['panelMass'])

        aircraftParameters['totalMass'] = aircraftParameters['rotorMass'] + aircraftParameters['panelMass'] + \
                                      aircraftParameters['batteryMass'] + aircraftParameters['wingMass'] + \
                                      aircraftParameters['tailMass'] + aircraftParameters['bodyMass'] + \
                                      const['payloadMass']
        # Payload-Range
        PayloadRange(aircraftParameters['rotorRadius'], aircraftParameters['rotorMass'], Mass_design,
                     const['ultimateLoad'], aircraftParameters['AR'], aircraftParameters['wingbraced'],
                     const['cruiseSpeed'], const['batteryEnergyDensity'], const['batteryPowerDensity'],
                     const['batteryEnergyDensity'], aircraftParameters['panelMass'], const['payloadMass'],
                     const['designRange'] / 1000, Mass_design - aircraftParameters['totalMass'])

        # Climb Performance
        ROC_cruise = AircraftClimbPerf(aircraftParameters['batteryMass'], const['batteryPowerDensity'], Mass_design,
                                       aircraftParameters['rotorRadius'], const['cruiseSpeed'])

        ROC_rotor = RotorClimbPerf(Mass_design, aircraftParameters['rotorRadius'],
                                   aircraftParameters['totalRotors'])

        print(f'Total aircraft weight = {aircraftParameters["totalMass"]}kg\n')
        print('-' * 50 + '\n')

        # Store results
        massArr = np.array([[aircraftParameters['totalMass'], aircraftParameters['rotorMass'],
                             aircraftParameters['wingMass'], aircraftParameters['tailMass'],
                             aircraftParameters['bodyMass'], aircraftParameters['batteryMass'],
                             aircraftParameters['panelMass']]]).T

        dimArr = np.array([[aircraftParameters['rotorRadius'], aircraftParameters['wingspan'],
                            aircraftParameters['chord']]]).T

        if iterate:
            diff = abs((aircraftParameters['totalMass'] - designMass) / designMass)
        else:
            diff = 0

    return massArr, dimArr, aircraftParameters



if __name__ == '__main__':
    global const, aircraftParameters

    # Sizing of Rotor:
    Mass_design = const['maxMass']  # kg
    v_arr = np.linspace(const['cruiseSpeed'], 130, 10)
    n_iterations = []

    con = const.copy()
    airc = aircraftParameters.copy()


    m, dim, acParams = design(iterate=True)
