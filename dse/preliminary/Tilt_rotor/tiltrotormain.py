from PerformanceAnalysis import PayloadRange, AircraftClimbPerf, RotorClimbPerf
from AircraftEstimating import Class2Weight, DragEstimation
from RotorEngineSizing import RadiusMassElementMomentum
from constants import const, aircraftParameters
from power_sizing import size_power_subsystem
from cruise_sizing import area_and_thrust
import xarray as xr
import numpy as np


if __name__ == '__main__':
    global const, aircraftParameters

    # Sizing of Rotor:
    Mass_design = const['maxMass']  # kg
    v_arr = np.linspace(const['cruiseSpeed'], 130, 10)
    n_iterations = []

    con = const.copy()
    airc = aircraftParameters.copy()

    # Set up Calculations
    for v in v_arr:
        const = con.copy()
        aircraftParameters = airc.copy()
        const['cruiseSpeed'] = v
        tipSpeed = -0.88*const['cruiseSpeed'] + 268.87
        diff = 100
        i = 0
        try:
            while diff > 0.01:
                Mass_thrust = Mass_design * (1 + const['margin'])  # kg
                aircraftParameters['rotorRadius'], takeOffThrustPerEngine, aircraftParameters['horsepowerPerEngine'], \
                aircraftParameters['totalPower'], aircraftParameters['rotorMass'] = \
                    RadiusMassElementMomentum(M=Mass_thrust, N_rotors=aircraftParameters['totalRotors'],
                                              N_blades=aircraftParameters['bladesPerRotor'],
                                              coaxial=aircraftParameters['coaxial'], V_tip=tipSpeed, print_results=True)


                # Calculate wing area
                aircraftParameters['wingArea'], aircraftParameters['cruiseThrust'] = \
                    area_and_thrust(0, const['cl'], const['cd'], Mass_design, 0.5*const['airDensity']*const['cruiseSpeed']**2)

                aircraftParameters['wingspan'] = 3 * aircraftParameters['rotorRadius']
                aircraftParameters['chord'] = aircraftParameters['wingArea']/aircraftParameters['wingspan']
                aircraftParameters['AR'] = aircraftParameters['wingspan']/aircraftParameters['chord']

                if aircraftParameters['wingArea']/aircraftParameters['wingspan'] < aircraftParameters['ARmin']:
                    aircraftParameters['wingspan'] = np.sqrt(aircraftParameters['wingArea']*aircraftParameters['ARmin'])
                    aircraftParameters['chord'] = aircraftParameters['wingArea']/aircraftParameters['wingspan']


                print(f'Wing area = {aircraftParameters["wingArea"]}[m2]')
                print(f'Wingspan = {aircraftParameters["wingspan"]}[m]')
                print(f'Chord = {aircraftParameters["chord"]}[m]')
                print(f'AR = {aircraftParameters["AR"]}')



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
                             const['designRange']/1000, Mass_design-aircraftParameters['totalMass'])

                # Climb Performance
                ROC_cruise = AircraftClimbPerf(aircraftParameters['batteryMass'], const['batteryPowerDensity'], Mass_design,
                                               aircraftParameters['rotorRadius'], const['cruiseSpeed'])

                ROC_rotor = RotorClimbPerf(Mass_design, aircraftParameters['rotorRadius'], aircraftParameters['totalRotors'])

                diff = abs((aircraftParameters['totalMass']-Mass_design)/Mass_design)
                print(f'Total aircraft weight = {aircraftParameters["totalMass"]}kg\n')
                print('-'*50 + '\n')

                # Store results
                if i == 0:
                    massArr = np.array([[aircraftParameters['totalMass'], aircraftParameters['rotorMass'],
                                         aircraftParameters['wingMass'], aircraftParameters['tailMass'],
                                         aircraftParameters['bodyMass'], aircraftParameters['batteryMass'],
                                         aircraftParameters['panelMass']]]).T

                    dimArr = np.array([[aircraftParameters['rotorRadius'], aircraftParameters['wingspan'],
                                        aircraftParameters['chord']]]).T


                else:
                    massArr = np.hstack((massArr, np.array([[aircraftParameters['totalMass'],
                                                             aircraftParameters['rotorMass'],
                                                             aircraftParameters['wingMass'],
                                                             aircraftParameters['tailMass'],
                                                             aircraftParameters['bodyMass'],
                                                             aircraftParameters['batteryMass'],
                                                             aircraftParameters['panelMass']]]).T))

                    dimArr = np.hstack((dimArr, np.array([[aircraftParameters['rotorRadius'],
                                                          aircraftParameters['wingspan'],
                                                          aircraftParameters['chord']]]).T))


                Mass_design = aircraftParameters['totalMass']
                n_iterations.append(i)
                i += 1

            m = list(v_arr).index(v)
            if m == 0:
                massResults = np.empty((np.size(v_arr), np.shape(massArr)[0], np.shape(massArr)[1]))
                massResults[0] = massArr

                dimResults = np.empty((np.size(v_arr), np.shape(dimArr)[0], np.shape(dimArr)[1]))
                dimResults[0] = dimArr
            else:
                if np.shape(massArr)[1] > np.shape(massResults[0])[1]:
                    n = np.shape(massArr)[1] - np.shape(massResults[0])[1]
                    massResults = np.pad(massResults, [(0, 0), (0, 0), (0, n)])
                    dimResults = np.pad(dimResults, [(0, 0), (0, 0), (0, n)])
                elif np.shape(massArr)[1] < np.shape(massResults[0])[1]:
                    n = np.shape(massResults[0])[1] - np.shape(massArr)[1]
                    massArr = np.pad(massArr, [(0, 0), (0, n)])
                    dimArr = np.pad(dimArr, [(0, 0), (0, n)])

                massResults[m] = massArr
                dimResults[m] = dimArr
        except:
            break

    n_it = []
    for i in range(np.shape(massResults)[-1]):
        n_it.append(f'Iteration {i}')


    massDs = xr.Dataset(
        {
            'mass': (('cruiseSpeed', 'parameter', 'iteration'), massResults),
        },
        coords={
            'cruiseSpeed': v_arr,
            'parameter': ['totalMass', 'rotorMass', 'wingMass', 'tailMass', 'bodyMass', 'batteryMass', 'panelMass'],
            'iteration': n_it,
        },
    )

    dimDS = xr.Dataset(
        {
            'dimensions': (('cruiseSpeed', 'parameter', 'iteration'), dimResults),
        },
        coords={
            'cruiseSpeed': v_arr,
            'parameter': ['rotor radius', 'wingspan', 'chord'],
            'iteration': n_it,
        },
    )

    df0 = massDs.to_dataframe()
    df0.to_excel('mass_iteration_results.xlsx')

    df1 = dimDS.to_dataframe()
    df1.to_excel('dimension_iteration_results.xlsx')
