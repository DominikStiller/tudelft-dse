import unittest
from scipy.stats.stats import pearsonr
from numpy.testing import assert_allclose
import numpy as np

from dse.preliminary.Tilt_rotor.cruise_sizing import area
from dse.preliminary.Tilt_rotor.constants import const,aircraftParameters
from dse.preliminary.Tilt_rotor.RotorEngineSizing import RadiusMassElementMomentum
from dse.preliminary.Tilt_rotor.AircraftEstimating import DragEstimation
from dse.preliminary.Tilt_rotor.power_sizing import size_power_subsystem
from dse.preliminary.Tilt_rotor.AircraftEstimating import Class2Weight

class RadiusMassElementMomentum_Unittest(unittest.TestCase):
    def test_correlation_mass_thrust(self):
        take_off_mass_lst = list(np.arange(100,3100,100))
        Thrust_lst = []
        for item in take_off_mass_lst:
            takeOffTipSpeed = 0.92 * const['soundSpeed']
            Thrust = RadiusMassElementMomentum(item, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False)[1]
            Thrust_lst.append(Thrust)
        correlation = pearsonr(take_off_mass_lst,Thrust_lst)
        self.assertTrue(correlation[0] > 0)

    def test_correlation_mass_radius(self):
        take_off_mass_lst = list(np.arange(100,3100,100))
        Radius_lst = []
        for item in take_off_mass_lst:
            takeOffTipSpeed = 0.92 * const['soundSpeed']
            Radius = RadiusMassElementMomentum(item, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False)[0]
            Radius_lst.append(Radius)
        correlation = pearsonr(take_off_mass_lst,Radius_lst)
        self.assertTrue(correlation[0] > 0)

    def test_correlation_mass_Rotor_mass(self):
        take_off_mass_lst = list(np.arange(100,3100,100))
        Rotor_mass_lst = []
        for item in take_off_mass_lst:
            takeOffTipSpeed = 0.92 * const['soundSpeed']
            Rotor_mass = RadiusMassElementMomentum(item, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False)[4]
            Rotor_mass_lst.append(Rotor_mass)
        correlation = pearsonr(take_off_mass_lst,Rotor_mass_lst)
        self.assertTrue(correlation[0] > 0)

    def test_correlation_N_rotors_thrust(self):
        N_rotors_lst = list(np.arange(1,31,1))
        Thrust_lst = []
        for item in N_rotors_lst:
            takeOffTipSpeed = 0.92 * const['soundSpeed']
            Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
            Thrust = RadiusMassElementMomentum(M=Mass_thrust, N_rotors=item,
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False)[1]
            Thrust_lst.append(Thrust)
        correlation = pearsonr(N_rotors_lst,Thrust_lst)
        self.assertTrue(correlation[0] < 0)

    def test_correlation_density_radius(self):
        Density_lst = list(np.arange(0.01,0.02,0.001))
        Radius_lst = []
        for item in Density_lst:
            takeOffTipSpeed = 0.92 * const['soundSpeed']
            Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
            Radius = RadiusMassElementMomentum(M=Mass_thrust, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False,changeConstants=('airDensity',item))[0]
            Radius_lst.append(Radius)
        correlation = pearsonr(Density_lst,Radius_lst)
        self.assertTrue(correlation[0] < 0)

    def test_correlation_gravity_thrust(self):
        Gravity_lst = list(np.arange(3,4,0.01))
        Thrust_lst = []
        for item in Gravity_lst:
            takeOffTipSpeed = 0.92 * const['soundSpeed']
            Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
            Thrust = RadiusMassElementMomentum(M=Mass_thrust, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False,changeConstants=('gravityMars',item))[1]
            Thrust_lst.append(Thrust)
        correlation = pearsonr(Gravity_lst,Thrust_lst)
        self.assertTrue(correlation[0] > 0)

    def test_zero_mass(self):
        mass = 0
        takeOffTipSpeed = 0.92 * const['soundSpeed']
        Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
        Radius, Thrust = RadiusMassElementMomentum(mass, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False, ExtremeValue=True)[:2]
        self.assertTrue(Radius == 1 and Thrust == 0)

    def test_zero_N_rotors(self):
        N_rotors = 0
        takeOffTipSpeed = 0.92 * const['soundSpeed']
        Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
        Return = RadiusMassElementMomentum(M=Mass_thrust, N_rotors=N_rotors,
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False, ExtremeValue=True)[0]
        self.assertTrue(Return == "N_rotors has to be greater than zero.")

    def test_zero_N_blades(self):
        N_blades = 0
        takeOffTipSpeed = 0.92 * const['soundSpeed']
        Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
        Return = RadiusMassElementMomentum(M=Mass_thrust, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=N_blades,
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False, ExtremeValue=True)[0]
        self.assertTrue(Return == "N_blades has to be greater than zero.")

class area_unittest(unittest.TestCase):
    def test_correlation_cl_area(self):
        cl_lst = list(np.arange(0.5,3,0.1))
        area_lst = []
        aircraftParameters['totalMass'] = const['maxMass']
        Mass_design = aircraftParameters['totalMass']
        for item in cl_lst:
            area(item, Mass_design, 0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2,Print=False)
            area_lst.append(aircraftParameters["wingArea"])
        correlation = pearsonr(cl_lst,area_lst)
        self.assertTrue(correlation[0] < 0)

    def test_correlation_mass_area(self):
        mass_lst = list(np.arange(100,3100,100))
        area_lst = []
        aircraftParameters['totalMass'] = const['maxMass']
        Mass_design = aircraftParameters['totalMass']
        for item in mass_lst:
            area(const['cl'], item, 0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2,Print=False)
            area_lst.append(aircraftParameters["wingArea"])
        correlation = pearsonr(mass_lst,area_lst)
        self.assertTrue(correlation[0] > 0)

    def test_correlation_q_area(self):
        q_lst = list(np.arange(10,300,10))
        area_lst = []
        aircraftParameters['totalMass'] = const['maxMass']
        Mass_design = aircraftParameters['totalMass']
        for item in q_lst:
            area(const['cl'], Mass_design, item, Print=False)
            area_lst.append(aircraftParameters["wingArea"])
        correlation = pearsonr(q_lst,area_lst)
        self.assertTrue(correlation[0] < 0)

    def test_correlation_gravity_area(self):
        gravity_lst = list(np.arange(3,4,0.01))
        area_lst = []
        aircraftParameters['totalMass'] = const['maxMass']
        Mass_design = aircraftParameters['totalMass']
        for item in gravity_lst:
            area(const['cl'], Mass_design, 0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2, Print=False, changeConstants=('gravityMars',item))
            area_lst.append(aircraftParameters["wingArea"])
        correlation = pearsonr(gravity_lst,area_lst)
        self.assertTrue(correlation[0] > 0)

    def test_zero_cl(self):
        cl = 0
        aircraftParameters['totalMass'] = const['maxMass']
        Mass_design = aircraftParameters['totalMass']
        area(cl, Mass_design, 0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2, Print=False)
        self.assertTrue(aircraftParameters["wingArea"] == np.infty)

    def test_zero_mass(self):
        mass = 0
        aircraftParameters['totalMass'] = const['maxMass']
        Mass_design = aircraftParameters['totalMass']
        area(const['cl'], mass, 0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2, Print=False)
        self.assertTrue(aircraftParameters["wingArea"] == 0.0)

class DragEstimation_Unittest(unittest.TestCase):
    def test_correlation_fuselage_size_cd0(self):
        lf_size_lst = list(np.arange(1,2.5,0.1))
        bf_size_lst = list(np.arange(1,2.5,0.1))
        cd01_lst = []
        cd02_lst = []
        for item in lf_size_lst:
            bf = 1
            cd0 = DragEstimation(aircraftParameters['wingArea'], const['cruiseSpeed'], const['visc_cr'],changeDimensions=(item,bf))
            cd01_lst.append(cd0)
        for item in bf_size_lst:
            lf = 1.78 + aircraftParameters['chord']
            cd0 = DragEstimation(aircraftParameters['wingArea'], const['cruiseSpeed'], const['visc_cr'],changeDimensions=(lf,item))
            cd02_lst.append(cd0)
        correlation_lf = pearsonr(lf_size_lst,cd01_lst)
        correlation_bf = pearsonr(bf_size_lst,cd02_lst)
        self.assertTrue(correlation_lf[0] > 0 and correlation_bf[0] > 0)

class size_power_subsystem_Unittest(unittest.TestCase):
    def test_correlation_T0thrust_battery_mass(self):
        T0thrust_lst = list(np.arange(2000,3000,100))
        battery_mass_lst = []

        takeOffTipSpeed = 0.92 * const['soundSpeed']
        Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
        aircraftParameters['rotorRadius'], takeOffThrustPerEngine, aircraftParameters['horsepowerPerEngine'], \
        aircraftParameters['totalPower'], aircraftParameters['rotorMass'] = \
            RadiusMassElementMomentum(M=Mass_thrust, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False)

        aircraftParameters['totalMass'] = const['maxMass']
        Mass_design = aircraftParameters['totalMass']
        area(const['cl'], Mass_design, 0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2, Print=False)

        cd0 = DragEstimation(aircraftParameters['wingArea'], const['cruiseSpeed'], const['visc_cr'])
        aircraftParameters['cruiseThrust'] = (cd0 + const['cl'] ** 2 / (
                np.pi * aircraftParameters['AR'] * const['oswald'])) * \
                                             0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2 * aircraftParameters[
                                                 'wingArea']
        const['airDensity'] = 0.01

        for item in T0thrust_lst:
            battery_mass = size_power_subsystem(aircraftParameters['rotorRadius'], item,
                                 aircraftParameters['cruiseThrust'],
                                 const['designRange'] / const['cruiseSpeed'],
                                 const['takeOffTime'], aircraftParameters['wingArea'], plot=False, Print=False)[0]
            battery_mass_lst.append(battery_mass)
        correlation = pearsonr(T0thrust_lst,battery_mass_lst)
        self.assertTrue(correlation[0] > 0)

    def test_correlation_time_battery_mass(self):
        cruise_time_lst = list(np.arange(1000,10000,200))
        takeoff_time_lst = list(np.arange(100,1000,100))
        battery_mass1_lst = []
        battery_mass2_lst = []

        takeOffTipSpeed = 0.92 * const['soundSpeed']
        Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
        aircraftParameters['rotorRadius'], takeOffThrustPerEngine, aircraftParameters['horsepowerPerEngine'], \
        aircraftParameters['totalPower'], aircraftParameters['rotorMass'] = \
            RadiusMassElementMomentum(M=Mass_thrust, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False)

        aircraftParameters['totalMass'] = const['maxMass']
        Mass_design = aircraftParameters['totalMass']
        area(const['cl'], Mass_design, 0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2, Print=False)

        cd0 = DragEstimation(aircraftParameters['wingArea'], const['cruiseSpeed'], const['visc_cr'])
        aircraftParameters['cruiseThrust'] = (cd0 + const['cl'] ** 2 / (
                    np.pi * aircraftParameters['AR'] * const['oswald'])) * \
                                             0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2 * aircraftParameters[
                                                 'wingArea']
        const['airDensity'] = 0.01

        for item in cruise_time_lst:
            battery_mass = size_power_subsystem(aircraftParameters['rotorRadius'], takeOffThrustPerEngine,
                                 aircraftParameters['cruiseThrust'], item,
                                 const['takeOffTime'], aircraftParameters['wingArea'], plot=False, Print=False)[0]
            battery_mass1_lst.append(battery_mass)
        for item in takeoff_time_lst:
            battery_mass = size_power_subsystem(aircraftParameters['rotorRadius'], takeOffThrustPerEngine,
                                 aircraftParameters['cruiseThrust'],
                                 const['designRange'] / const['cruiseSpeed'],
                                 item, aircraftParameters['wingArea'], plot=False, Print=False)[0]
            battery_mass2_lst.append(battery_mass)
        correlation1 = pearsonr(cruise_time_lst,battery_mass1_lst)
        correlation2 = pearsonr(takeoff_time_lst,battery_mass2_lst)
        self.assertTrue(correlation1[0] > 0 and correlation2[0] > 0)

    def test_power_surplus_high_latitudes(self):
        takeOffTipSpeed = 0.92 * const['soundSpeed']
        Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
        aircraftParameters['rotorRadius'], takeOffThrustPerEngine, aircraftParameters['horsepowerPerEngine'], \
        aircraftParameters['totalPower'], aircraftParameters['rotorMass'] = \
            RadiusMassElementMomentum(M=Mass_thrust, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False)

        aircraftParameters['totalMass'] = const['maxMass']
        Mass_design = aircraftParameters['totalMass']
        area(const['cl'], Mass_design, 0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2, Print=False)

        cd0 = DragEstimation(aircraftParameters['wingArea'], const['cruiseSpeed'], const['visc_cr'])
        aircraftParameters['cruiseThrust'] = (cd0 + const['cl'] ** 2 / (
                np.pi * aircraftParameters['AR'] * const['oswald'])) * \
                                             0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2 * aircraftParameters[
                                                 'wingArea']
        const['airDensity'] = 0.01

        powerSurplus = size_power_subsystem(aircraftParameters['rotorRadius'], takeOffThrustPerEngine,
                                 aircraftParameters['cruiseThrust'],
                                 const['designRange'] / const['cruiseSpeed'],
                                 const['takeOffTime'], aircraftParameters['wingArea'], plot=False, Print=False)[2]
        self.assertTrue(np.average(powerSurplus[60]) < np.average(powerSurplus[120]))

class Class2weight_Unittest(unittest.TestCase):
    def test_correlation_rotor_radius_tail_area(self):
        rotor_radius_lst = list(np.arange(1,25,1))
        tail_area_lst = []

        takeOffTipSpeed = 0.92 * const['soundSpeed']
        Mass_thrust = const['maxMass'] * (1 + const['margin'])  # kg
        aircraftParameters['rotorRadius'], takeOffThrustPerEngine, aircraftParameters['horsepowerPerEngine'], \
        aircraftParameters['totalPower'], aircraftParameters['rotorMass'] = \
            RadiusMassElementMomentum(M=Mass_thrust, N_rotors=aircraftParameters['totalRotors'],
                                      N_blades=aircraftParameters['bladesPerRotor'],
                                      coaxial=aircraftParameters['coaxial'], V_tip=takeOffTipSpeed, print_results=False)
        aircraftParameters['totalMass'] = const['maxMass']
        Mass_design = aircraftParameters['totalMass']
        area(const['cl'], Mass_design, 0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2, Print=False)
        cd0 = DragEstimation(aircraftParameters['wingArea'], const['cruiseSpeed'], const['visc_cr'])
        aircraftParameters['cruiseThrust'] = (cd0 + const['cl'] ** 2 / (
                    np.pi * aircraftParameters['AR'] * const['oswald'])) * \
                                             0.5 * const['airDensity'] * const['cruiseSpeed'] ** 2 * aircraftParameters[
                                                 'wingArea']
        const['airDensity'] = 0.01
        aircraftParameters['batteryMass'], aircraftParameters['panelMass'], powerSurplus = \
            size_power_subsystem(aircraftParameters['rotorRadius'], takeOffThrustPerEngine,
                                 aircraftParameters['cruiseThrust'],
                                 const['designRange'] / const['cruiseSpeed'],
                                 const['takeOffTime'], aircraftParameters['wingArea'], plot=False, Print=False)

        for item in rotor_radius_lst:
            tail_area = Class2Weight(item, aircraftParameters['rotorMass'], Mass_design,
                         const['ultimateLoad'], aircraftParameters['AR'], aircraftParameters['wingbraced'],
                         const['cruiseSpeed'], const['takeoffBatteryEnergyDensity'], const['takeoffBatteryEnergyDensity'],
                         const['takeoffBatteryPowerDensity'], const['payloadMass'], aircraftParameters['panelMass'], print_results=False)[2]
            tail_area_lst.append(tail_area)

        correlation = pearsonr(rotor_radius_lst,tail_area_lst)
        self.assertTrue(correlation[0] < 0)


if __name__ == '__main__':
    unittest.main()