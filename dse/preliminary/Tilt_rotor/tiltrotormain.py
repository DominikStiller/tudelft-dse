from PerformanceAnalysis import PayloadRange, AircraftClimbPerf, RotorClimbPerf
from AircraftEstimating import Class2Weight, DragEstimation
from RotorEngineSizing import RadiusMassElementMomentum
from structures import calculate_cg, assembly_volume
from constants import const, aircraftParameters
from power_sizing import size_power_subsystem
from cruise_sizing import area_and_thrust
from matplotlib import pyplot as plt
import numpy as np


def design(iterate=False):
    diff = 100
    Range = 1000
    aircraftParameters['totalMass'] = const['maxMass']
    while diff > 0.01 and Range >= 1000:
        Mass_design = aircraftParameters['totalMass']
        tipSpeed = -0.88 * const['cruiseSpeed'] + 268.87
        Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
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
            print(f'Wing span is too large for the required AR')
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
        print(f'Total Drag on the aircraft: {aircraftParameters["cruiseThrust"]}[N]')
        aircraftParameters['batteryMass'], aircraftParameters['panelMass'], powerSurplus = \
            size_power_subsystem(aircraftParameters['rotorRadius'], takeOffThrustPerEngine,
                                 aircraftParameters['cruiseThrust'],
                                 const['designRange'] / const['cruiseSpeed'] + const['takeOffTime'],
                                 const['takeOffTime'], aircraftParameters['wingArea'], plot=False)

        # Calculate weights
        Range, aircraftParameters['wingMass'], aircraftParameters['tailMass'], aircraftParameters['bodyMass'] = \
            Class2Weight(aircraftParameters['rotorRadius'], aircraftParameters['rotorMass'], Mass_design,
                         const['ultimateLoad'], aircraftParameters['AR'], aircraftParameters['wingbraced'],
                         const['cruiseSpeed'], const['batteryEnergyDensity'], const['batteryPowerDensity'],
                         const['batteryEnergyDensity'], const['payloadMass'], aircraftParameters['panelMass'])

        print(f'Maximum Range available: {Range}[km]')

        aircraftParameters['totalMass'] = aircraftParameters['rotorMass'] + aircraftParameters['panelMass'] + \
                                      aircraftParameters['batteryMass'] + aircraftParameters['wingMass'] + \
                                      aircraftParameters['tailMass'] + aircraftParameters['bodyMass'] + \
                                      const['payloadMass']

        if aircraftParameters['totalMass'] > Mass_design:
            print('Total mass was higher than the design mass')
            return 0, 0, 0

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

        aircraftParameters['aircraftVolume'] = assembly_volume()
        print(f'Total aircraft weight = {aircraftParameters["totalMass"]}kg\n')
        print(f'Cg location is {1.78 + aircraftParameters["chord"] / 2 - calculate_cg()} in front of the rotors thrust')
        print(f'Aircraft volume = {aircraftParameters["aircraftVolume"]}')
        print('-' * 50 + '\n')

        # Store results
        massArr = np.array([[aircraftParameters['totalMass'], aircraftParameters['rotorMass'],
                             aircraftParameters['wingMass'], aircraftParameters['tailMass'],
                             aircraftParameters['bodyMass'], aircraftParameters['batteryMass'],
                             aircraftParameters['panelMass']]]).T

        dimArr = np.array([[aircraftParameters['rotorRadius'], aircraftParameters['wingspan'],
                            aircraftParameters['chord']]]).T


        if iterate:
            diff = abs((aircraftParameters['totalMass'] - Mass_design) / Mass_design)
        else:
            diff = 0

    return massArr, dimArr, aircraftParameters


if __name__ == '__main__':
    global const, aircraftParameters

    # Sizing of Rotor:
    con = const.copy()
    airc = aircraftParameters.copy()


    m, dim, acParams = design(iterate=True)

    plt.legend(loc='best')
    plt.savefig('../Figures/Payload-Range.pdf')
    plt.show()
