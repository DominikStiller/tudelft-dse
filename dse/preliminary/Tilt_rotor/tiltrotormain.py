from PerformanceAnalysis import PayloadRange, AircraftClimbPerf, RotorClimbPerf
from AircraftEstimating import Class2Weight, DragEstimation
from RotorEngineSizing import RadiusMassElementMomentum
from structures import calculate_cg, assembly_volume, max_rotor_loads
from constants import const, aircraftParameters
from power_sizing import size_power_subsystem
from cruise_sizing import area
from matplotlib import pyplot as plt
import numpy as np


def design(iterate=False):
    diff = 100
    Range = 1000
    aircraftParameters['totalMass'] = const['maxMass']
    while diff > 0.005 and Range >= 1000:
        const['airDensity'] = 0.01
        Mass_design = aircraftParameters['totalMass']

        takeOffTipSpeed = 0.92 * const['soundSpeed']
        Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
        aircraftParameters['rotorRadius'], takeOffThrustPerEngine, aircraftParameters['horsepowerPerEngine'], \
        aircraftParameters['totalPower'], aircraftParameters['rotorMass'] = \
            RadiusMassElementMomentum(M=Mass_thrust, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=True)

        if aircraftParameters['rotorRadius'] == "N_rotors has to be greater than zero.":
            raise ValueError("N_rotors has to be greater than zero.")
        elif aircraftParameters['rotorRadius'] == "N_blades has to be greater than zero.":
            raise ValueError("N_blades has to be greater than zero.")

        # Calculate wing area
        area(const['cl'], Mass_design, 0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2)

        cd0 = DragEstimation(aircraftParameters['wingArea'], const['cruiseSpeed'], const['visc_cr'])
        aircraftParameters['cruiseThrust'] = (cd0 + const['cl']**2/(np.pi*aircraftParameters['AR']*const['oswald'])) * \
            0.5 * const['airDensity'] * const['cruiseSpeed']**2 * aircraftParameters['wingArea']
        print(f'Cruise drag = {aircraftParameters["cruiseThrust"]}')

        const['airDensity'] = 0.01

        # Size batteries for cruise and solar panels
        print(f'Total Drag on the aircraft: {aircraftParameters["cruiseThrust"]}[N]')
        aircraftParameters['batteryMass'], aircraftParameters['panelMass'], powerSurplus = \
            size_power_subsystem(aircraftParameters['rotorRadius'], takeOffThrustPerEngine,
                                 aircraftParameters['cruiseThrust'],
                                 const['designRange'] / const['cruiseSpeed'],
                                 const['takeOffTime'], aircraftParameters['wingArea'], plot=False)

        # Calculate weights
        Range, aircraftParameters['wingMass'], aircraftParameters['tailMass'], aircraftParameters['bodyMass'] = \
            Class2Weight(aircraftParameters['rotorRadius'], aircraftParameters['rotorMass'], Mass_design,
                         const['ultimateLoad'], aircraftParameters['AR'], aircraftParameters['wingbraced'],
                         const['cruiseSpeed'], const['takeoffBatteryEnergyDensity'], const['takeoffBatteryEnergyDensity'],
                         const['takeoffBatteryPowerDensity'], const['payloadMass'], aircraftParameters['panelMass'])

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
    s1223 = np.array([[1, 0.99838, 0.99417, 0.98825, 0.98075, 0.97111, 0.95884,0.94389,0.92639,0.90641,0.88406,0.85947,0.83277,
                  0.80412,0.77369,0.74166,0.70823,0.6736,0.63798,0.60158,0.56465,0.52744,0.49025,0.4534,0.41721,0.38193,
                  0.34777,0.31488,0.28347,0.2537,0.22541,0.19846,0.17286,0.14863,0.12591,0.10482,0.08545,0.06789,0.05223,
                  0.03855,0.02694,0.01755,0.01028,0.00495,0.00155,0.00005,0.00044,0.00264,0.00789,0.01718,0.03006,0.04627,
                  0.06561,0.08787,0.11282,0.1402,0.17006,0.20278,0.2384,0.27673,0.3175,0.36044,0.40519,0.45139,0.4986,0.54639,
                  0.59428,0.64176,0.68832,0.73344,0.7766,0.81729,0.855,0.88928,0.91966,0.94573,0.96693,0.98255,0.99268,0.99825,1],
             [0,0.00126,0.00494,0.01037,0.01646,0.0225,0.02853,0.03476,0.04116,0.04768,0.05427,0.06089,0.06749,0.07402,0.08044,
              0.08671,0.09277,0.09859,0.10412,0.10935,0.11425,0.11881,0.12303,0.12683,0.13011,0.13271,0.13447,0.13526,0.13505,
              0.13346,0.13037,0.12594,0.12026,0.11355,0.10598,0.0977,0.08879,0.0794,0.06965,0.05968,0.04966,0.03961,0.02954,0.01969,
              0.01033,0.00178,-0.00561,-0.0112,-0.01427,-0.0155,-0.01584,-0.01532,-0.01404,-0.01202,-0.00925,-0.00563,-0.00075,0.00535,
              0.01213,0.01928,0.02652,0.03358,0.04021,0.04618,0.05129,0.05534,0.0582,0.05976,0.05994,0.05872,0.05612,0.05219,0.04706,
              0.04088,0.03387,0.02624,0.01822,0.0106,0.00468,0.00115,0]])
    max_rotor_loads(aircraftParameters['rotorRadius']/20 *s1223)

    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          ncol=3, fancybox=True, shadow=True)
    plt.savefig('../Figures/Payload-Range.pdf')
    plt.show()
