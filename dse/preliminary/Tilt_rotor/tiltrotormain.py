from PerformanceAnalysis import PayloadRange, AircraftClimbPerf, RotorClimbPerf
from AircraftEstimating import Class2Weight, DragEstimation
from RotorEngineSizing import RadiusMassElementMomentum
from constants import const, aircraftParameters
from power_sizing import size_power_subsystem
from cruise_sizing import area_and_thrust
import xarray as xr
import numpy as np


def design(constants, aircraftParams, iterate=False):
    diff = 100
    while diff > 0.01:
        designMass = aircraftParams['totalMass']
        tipSpeed = -0.88 * constants['cruiseSpeed'] + 268.87
        Mass_thrust = Mass_design * (1 + constants['margin'])  # kg
        aircraftParams['rotorRadius'], takeOffThrustPerEngine, aircraftParams['horsepowerPerEngine'], \
        aircraftParams['totalPower'], aircraftParams['rotorMass'] = \
            RadiusMassElementMomentum(M=Mass_thrust, N_rotors=aircraftParams['totalRotors'],
                                      N_blades=aircraftParams['bladesPerRotor'],
                                      coaxial=aircraftParams['coaxial'], V_tip=tipSpeed, print_results=True)

        # Calculate wing area
        aircraftParams['wingArea'], aircraftParams['cruiseThrust'] = \
            area_and_thrust(0, constants['cl'], const['cd'], Mass_design,
                            0.5 * constants['airDensity'] * constants['cruiseSpeed'] ** 2)

        aircraftParams['wingspan'] = max(3 * aircraftParams['rotorRadius'],
                                         np.sqrt(aircraftParams['wingArea'] * aircraftParams['AR']))
        aircraftParams['chord'] = aircraftParams['wingspan'] / aircraftParams['AR']
        aircraftParams['wingArea'] = aircraftParams['wingspan'] * aircraftParams['chord']

        print(f'Wing area = {aircraftParams["wingArea"]}[m2]')
        print(f'Wingspan = {aircraftParams["wingspan"]}[m]')
        print(f'Chord = {aircraftParams["chord"]}[m]')
        print(f'AR = {aircraftParams["wingspan"] / aircraftParams["chord"]}')

        # Size batteries for cruise and solar panels
        aircraftParams['cruiseThrust'] += DragEstimation(R=aircraftParams['rotorRadius'],
                                                         Swing=aircraftParams['wingArea'], t2c=constants['t/c'],
                                                         Vcr=constants['cruiseSpeed'], visc_cr=constants['visc_cr'],
                                                         AR=aircraftParams['AR'])

        aircraftParams['batteryMass'], aircraftParams['panelMass'], powerSurplus = \
            size_power_subsystem(aircraftParams['rotorRadius'], takeOffThrustPerEngine,
                                 aircraftParams['cruiseThrust'],
                                 constants['designRange'] / constants['cruiseSpeed'] + constants['takeOffTime'],
                                 constants['takeOffTime'], aircraftParams['wingArea'], plot=True)

        # Calculate weights
        Range, aircraftParams['wingMass'], aircraftParams['tailMass'], aircraftParams['bodyMass'] = \
            Class2Weight(aircraftParams['rotorRadius'], aircraftParams['rotorMass'], constants['maxMass'],
                         constants['ultimateLoad'], aircraftParams['AR'], aircraftParams['wingbraced'],
                         constants['cruiseSpeed'], constants['batteryEnergyDensity'], constants['batteryPowerDensity'],
                         constants['batteryEnergyDensity'], constants['payloadMass'], aircraftParams['panelMass'])

        aircraftParams['totalMass'] = aircraftParams['rotorMass'] + aircraftParams['panelMass'] + \
                                      aircraftParams['batteryMass'] + aircraftParams['wingMass'] + \
                                      aircraftParams['tailMass'] + aircraftParams['bodyMass'] + \
                                      constants['payloadMass']
        # Payload-Range
        PayloadRange(aircraftParams['rotorRadius'], aircraftParams['rotorMass'], Mass_design,
                     constants['ultimateLoad'], aircraftParams['AR'], aircraftParams['wingbraced'],
                     constants['cruiseSpeed'], constants['batteryEnergyDensity'], constants['batteryPowerDensity'],
                     constants['batteryEnergyDensity'], aircraftParams['panelMass'], constants['payloadMass'],
                     constants['designRange'] / 1000, Mass_design - aircraftParams['totalMass'])

        # Climb Performance
        ROC_cruise = AircraftClimbPerf(aircraftParams['batteryMass'], constants['batteryPowerDensity'], Mass_design,
                                       aircraftParams['rotorRadius'], constants['cruiseSpeed'])

        ROC_rotor = RotorClimbPerf(Mass_design, aircraftParams['rotorRadius'],
                                   aircraftParams['totalRotors'])

        print(f'Total aircraft weight = {aircraftParams["totalMass"]}kg\n')
        print('-' * 50 + '\n')

        # Store results
        massArr = np.array([[aircraftParams['totalMass'], aircraftParams['rotorMass'],
                             aircraftParams['wingMass'], aircraftParams['tailMass'],
                             aircraftParams['bodyMass'], aircraftParams['batteryMass'],
                             aircraftParams['panelMass']]]).T

        dimArr = np.array([[aircraftParams['rotorRadius'], aircraftParams['wingspan'],
                            aircraftParams['chord']]]).T

        if iterate:
            diff = abs((aircraftParams['totalMass'] - designMass) / designMass)
        else:
            diff = 0

    return massArr, dimArr, aircraftParams



if __name__ == '__main__':
    global const, aircraftParameters

    # Sizing of Rotor:
    Mass_design = const['maxMass']  # kg
    v_arr = np.linspace(const['cruiseSpeed'], 130, 10)
    n_iterations = []

    con = const.copy()
    airc = aircraftParameters.copy()


    m, dim, acParams = design(con, airc, iterate=True)
