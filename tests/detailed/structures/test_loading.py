from unittest import TestCase
import numpy as np
import os
from dse.detailed.Structures.StructureClasses import *
from dse.detailed.Structures.material_properties import materials


class Test(TestCase):
    def test_xflr_forces(self):
        from dse.detailed.Structures.StructureClasses import xflr_forces

        if os.getcwd().split("\\")[-1] == "structures":
            os.chdir("..\\..\\..\\dse\\detailed\\Structures")
        elif os.getcwd().split("\\")[-1] != "Structures":
            os.chdir("..\\dse\\detailed\\Structures")

        force = xflr_forces(filename="Test_xflr5_file.csv", q=0.5 * 0.01 * 112 ** 2, b=16.8)

        if np.all(force.fy == 0):
            if np.all(force.fx < 0):
                if np.all(force.fz > 0):
                    if np.shape(force.F) == np.shape(force.application):
                        assert True
                    else:
                        assert np.shape(force.F) == np.shape(force.application)
                else:
                    assert np.all(force.fz < 0)
            else:
                assert np.all(force.fx > 0)
        else:
            assert np.all(force.fy == 0)

    def test_xflr_forces_file_type(self):
        from dse.detailed.Structures.StructureClasses import xflr_forces

        try:
            Forces = xflr_forces(5, q=0.5 * 0.01 * 112 ** 2, b=16.8)
            a = False
        except TypeError:
            a = True
        self.assertTrue(a)

    def test_xflr_forces_q_type(self):
        from dse.detailed.Structures.StructureClasses import xflr_forces

        if os.getcwd().split("\\")[-1] == "structures":
            os.chdir("..\\..\\..\\dse\\detailed\\Structures")
        elif os.getcwd().split("\\")[-1] != "Structures":
            os.chdir("..\\dse\\detailed\\Structures")

        try:
            Forces = xflr_forces("Test_xflr5_file.csv", q="0.5 * 0.01 * 112 ** 2", b=16.8)
            a = False
        except TypeError:
            a = True
        self.assertTrue(a)

    def test_xflr_forces_b_type(self):
        from dse.detailed.Structures.StructureClasses import xflr_forces

        if os.getcwd().split("\\")[-1] == "structures":
            os.chdir("..\\..\\..\\dse\\detailed\\Structures")
        elif os.getcwd().split("\\")[-1] != "Structures":
            os.chdir("..\\dse\\detailed\\Structures")

        try:
            Forces = xflr_forces(filename="Test_xflr5_file.csv", q=0.5 * 0.01 * 112 ** 2, b="16.8")
            a = False
        except TypeError:
            a = True
        self.assertTrue(a)




class TestForce(TestCase):
    def test_concatenate_two_forces(self):
        from dse.detailed.Structures.StructureClasses import Force

        dummy = Force(
            magnitude=np.array([[3, 3, 3]]).T, point_of_application=np.array([[1, 0, 0]]).T
        )
        dummy_2 = Force(
            magnitude=np.array([[-2, 0, 1]]).T, point_of_application=np.array([[2, 0, 0]]).T
        )
        dummy = dummy.concatenate_forces(dummy_2)

        assert np.shape(dummy.F) == (3, 2)

    def test_concatenate_functions_list(self):
        from dse.detailed.Structures.StructureClasses import Force

        dummy = Force(
            magnitude=np.array([[3, 3, 3]]).T, point_of_application=np.array([[1, 0, 0]]).T
        )
        dummy_2 = Force(
            magnitude=np.array([[-2, 0, 1]]).T, point_of_application=np.array([[2, 0, 0]]).T
        )
        dummy_3 = Force(
            magnitude=np.array([[1, 0, 0]]).T, point_of_application=np.array([[3, 0, 0]]).T
        )

        dummy = dummy.concatenate_forces([dummy_2, dummy_3])

        assert np.shape(dummy.F) == (3, 3)

    def test_separate_forces_number(self):
        from dse.detailed.Structures.StructureClasses import Force

        mag = np.array([[1, 5], [0, 0], [0, 0]])
        app = np.array([[1, 3], [0, 0], [0, 0]])
        dummy = Force(magnitude=mag, point_of_application=app)
        lst = dummy.separate_forces()

        if len(lst) == np.shape(mag)[1]:
            for i in range(len(lst)):
                if np.shape(lst[i].F) != (3, 1):
                    assert np.shape(lst[i].F) == (3, 1)
            assert True
        else:
            assert len(lst) == np.shape(mag)[1]


class TestBeam(TestCase):
    def test_add_loading_point_load(self):
        from dse.detailed.Structures.StructureClasses import Beam, Force

        x = 0.1
        y = 2
        z = 0.05
        wing = Beam(x, y, z, "square", material='Al/Si',
                    fixing_points=np.array([[x/2], [z/2]]) * np.ones(100))

        F = 1
        loc = -1.5

        point_load = Force(
            magnitude=np.array([[0, 0, F]]).T, point_of_application=np.array([[x/2, loc, z/2]]).T
        )

        wing.add_loading(point_load)

        y_index = (np.abs(wing.y - loc)).argmin()
        if np.all(wing.f_loading[y_index:] == np.array([[0, 0, F]]).T):
            if (np.abs((wing.m_loading[y_index] - wing.m_loading[-1]) / loc)[0] - F) < 0.01:
                assert True
            else:
                assert (np.abs((wing.m_loading[y_index] - wing.m_loading[-1]) / loc)[0] - F) < 0.01
        else:
            assert np.all(wing.f_loading[y_index:] == np.array([[0, 1, 0]]).T)

    def test_add_loading_distributed_load(self):
        from dse.detailed.Structures.StructureClasses import Beam, Force

        x = 0.1
        y = 2
        z = 0.05
        wing = Beam(x, y, z, "square", material='Al/Si',
                    fixing_points=np.array([[x/2], [z/2]]) * np.ones(100))

        F = np.array([[0, 0, 0], [0, 0, 0], [1, 1, 1]])

        app = np.array([[x/2, x/2, x/2], [-0.1, -1, -1.5], [z/2, z/2, z/2]])

        force = Force(magnitude=F, point_of_application=app)
        wing.add_loading(force)

        y_index_1 = (np.abs(wing.y - app[:, 0][1])).argmin()
        y_index_2 = (np.abs(wing.y - app[:, 1][1])).argmin()
        y_index_3 = (np.abs(wing.y - app[:, 2][1])).argmin()

        if (
                np.all(wing.f_loading[y_index_1:] == np.reshape(np.sum(F, 1), (3, 1)))
                and np.all(
            wing.f_loading[y_index_2:y_index_1] == np.reshape(np.sum(F[:, 1:], 1), (3, 1))
        )
                and np.all(wing.f_loading[y_index_3:y_index_2] == np.reshape(F[:, -1], (3, 1)))
                and np.all(wing.f_loading[:y_index_3] == 0)
        ):
            if (
                    np.abs(wing.m_loading[-1][0]) - 2.6 < 0.01
                    and np.abs(wing.m_loading[y_index_1][0]) - 2.3 < 0.01
                    and np.abs(wing.m_loading[y_index_2][0]) - 0.5 < 0.01
                    and np.all(wing.m_loading[:y_index_3] == 0)
            ):
                assert True
            else:
                assert np.abs(wing.m_loading[-1][0]) - 2.6 < 0.01, "Moment at the root is off"
                assert np.abs(wing.m_loading[y_index_1][0]) - 2.3 < 0.01, "Moment at point 1 is off"
                assert np.abs(wing.m_loading[y_index_2][0]) - 0.5 < 0.01, "Moment at point 2 is off"
                assert np.all(wing.m_loading[:y_index_3] == 0), "Moment at the tip is not 0"
        else:
            assert np.all(
                wing.f_loading[y_index_1:] == np.reshape(np.sum(F, 1), (3, 1))
            ), "Forces at root are off"
            assert np.all(
                wing.f_loading[y_index_2:y_index_1] == np.reshape(np.sum(F[:, 1:], 1), (3, 1))
            ), "Forces at middle are off"
            assert np.all(
                wing.f_loading[y_index_3:y_index_2] == np.reshape(F[:, -1], (3, 1))
            ), "Forces at semi-middle are off"
            assert np.all(wing.f_loading[:y_index_3] == 0), "Forces at the tip are off"

    def test_InternalLoads_Rhombus(self):
        ## Hand Calculations
        NAx_hand = 0.  # Neutral axis as obtained by hand calculations
        NAz_hand = 0.  # Neutral axis as obtained by hand calculations
        Ixx_hand = 0.0375  # [m4]
        Izz_hand = 0.0375  # [m4]
        Ixz_hand = 0.000   # [m4]
        sigma_hand = np.array(
            [
                [0.266666666666, 1.4666666666, 2.666666666, 1.2, -0.266666666666, -1.4666666666, -2.6666666, -1.2]
            ]
        ).T * 1e6
        b = 10
        x = np.array([[1, 0.5, 0, -0.5, -1, -0.5, 0, 0.5]])
        z = np.array([[0, 0.5, 1, 0.5, 0, -0.5, -1, -0.5]])
        l = np.linspace(-b, 0, 100)
        rombus = Beam(
            width=x.T * np.ones(100),
            height=z.T * np.ones(100),
            length=l,
            cross_section=np.vstack((x, z)) * np.ones((np.size(l), 2, 1)),
            material='Al/Si',
            fixing_points=np.array([[0], [0]]) * np.ones(100)
        )
        Applied_Load = Force(
            magnitude=np.array(
                [
                    [1000],
                    [0],
                    [-10000]
                ]
            ),
            point_of_application=np.array(
                [
                    [0],
                    [-b],
                    [0]
                ]
            )
        )
        rombus.add_loading(Applied_Load)
        # rombus.plot_internal_loading()

        x_booms_nr, z_booms_nr = np.split(np.reshape(rombus.section, (np.size(rombus.y), 2 * np.shape(rombus.x)[0])), 2,
                                          1)
        if np.all(x_booms_nr[:, 0] == x_booms_nr[:, -1]):
            x_booms_nr = x_booms_nr[:, :-1]
            z_booms_nr = x_booms_nr[:, :-1]

        x_booms_nr = x_booms_nr.T
        z_booms_nr = z_booms_nr.T

        NAx, NAz = rombus.NeutralAxis(boomArea_nr=np.ones((8, np.size(l)))*0.0125, x_booms_nr=x_booms_nr, z_booms_nr=z_booms_nr)
        Ixx, Izz, Ixz = rombus.MoI(boomArea_nr=np.ones((8, np.size(l)))*0.0125, x_booms_nr=x_booms_nr, z_booms_nr=z_booms_nr)
        sigma_nr = rombus.StressCalculations(boomArea_nr=np.ones((8, np.size(l)))*0.0125)

        assert abs((NAx[-1] - NAx_hand)) < 0.0001, 'The NAx is off'
        assert abs((NAz[-1] - NAz_hand)) < 0.0001, 'The NAz is off'
        assert abs(Ixx[-1] - Ixx_hand) < 0.0001, 'The Ixx is off'
        assert abs(Izz[-1] - Izz_hand) < 0.0001, 'The Izz is off'
        assert abs(Ixz[-1] - Ixz_hand) < 0.0001, 'The Ixz is off'
        assert np.all(abs(sigma_nr[:,-1] - sigma_hand.T)/sigma_hand.T < 0.01), 'The stresses are wrong'

    def test_InternalLoads_Parallelogram(self):
        ## Hand Calculations
        NAx_hand = 0.  # Neutral axis as obtained by hand calculations
        NAz_hand = 0.  # Neutral axis as obtained by hand calculations
        Ixx_hand = 0.0375  # [m4]
        Izz_hand = 0.075  # [m4]
        Ixz_hand = 0.0375  # [m4]
        sigma_hand = np.array(
            [
                [-2.4, 0.133333, 2.6666, 2.5333, 2.4, -1.3333, -2.666, -2.53333]
            ]
        ).T * 1e6
        b = 10
        x = np.array([[1, 1, 1, 0, -1, -1, -1, 0]])
        z = np.array([[0, 0.5, 1, 0.5, 0, -0.5, -1, -0.5]])
        l = np.linspace(-b, 0, 100)
        parallelogram = Beam(
            width=x.T * np.ones(np.size(l)),
            height=z.T * np.ones(np.size(l)),
            length=l,
            cross_section=np.vstack((x, z)) * np.ones((np.size(l), 2, 1)),
            material='Al/Si',
            fixing_points=np.array([[0], [0]]) * np.ones(np.size(l))
        )
        Applied_Load = Force(
            magnitude=np.array(
                [
                    [1000],
                    [0],
                    [-10000]
                ]
            ),
            point_of_application=np.array(
                [
                    [0],
                    [-b],
                    [0]
                ]
            )
        )
        parallelogram.add_loading(Applied_Load)

        x_booms_nr, z_booms_nr = np.split(np.reshape(parallelogram.section, (np.size(parallelogram.y), 2 * np.shape(parallelogram.x)[0])), 2,
                                          1)
        if np.all(x_booms_nr[:, 0] == x_booms_nr[:, -1]):
            x_booms_nr = x_booms_nr[:, :-1]
            z_booms_nr = x_booms_nr[:, :-1]

        x_booms_nr = x_booms_nr.T
        z_booms_nr = z_booms_nr.T

        NAx, NAz = parallelogram.NeutralAxis(boomArea_nr=np.ones((8, np.size(l)))*0.0125, x_booms_nr=x_booms_nr, z_booms_nr=z_booms_nr)
        Ixx, Izz, Ixz = parallelogram.MoI(boomArea_nr=np.ones((8, np.size(l)))*0.0125, x_booms_nr=x_booms_nr, z_booms_nr=z_booms_nr)
        sigma_nr = parallelogram.StressCalculations(boomArea_nr=np.ones((8, np.size(l)))*0.0125)

        assert abs((NAx[-1] - NAx_hand)) < 0.0001, 'The NAx is off'
        assert abs((NAz[-1] - NAz_hand)) < 0.0001, 'The NAz is off'
        assert abs(Ixx[-1] - Ixx_hand) < 0.0001, 'The Ixx is off'
        assert abs(Izz[-1] - Izz_hand) < 0.0001, 'The Izz is off'
        assert abs(Ixz[-1] - Ixz_hand) < 0.0001, 'The Ixz is off'
        assert np.all(abs(sigma_nr[:,-1] - sigma_hand.T)/sigma_hand.T < 0.01), 'The stresses are wrong'

    def test_neutral_axis(self):
        from dse.detailed.Structures.StructureClasses import Beam
        x = np.atleast_2d(np.hstack(
            (np.linspace(-1, 0, 25)[:-1],
             np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25)[:-1])
        ))
        z = np.atleast_2d(np.hstack(
            (np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25)[:-1],
             np.linspace(-1, 0, 25)[:-1]
             )
        ))
        rombus = Beam(
            width=x.T * np.ones(100),
            height=z.T * np.ones(100),
            length=5,
            cross_section=np.vstack((x, z)) * np.ones((100, 1, 1)),
            material='Al/Si',
            fixing_points=np.array([[0], [0]]) * np.ones(100)
        )
        A = np.ones((np.size(x), 100))

        x_booms_nr, z_booms_nr = np.split(np.reshape(rombus.section, (np.size(rombus.y), 2 * np.shape(rombus.x)[0])), 2,
                                          1)
        if np.all(x_booms_nr[:, 0] == x_booms_nr[:, -1]):
            x_booms_nr = x_booms_nr[:, :-1]
            z_booms_nr = x_booms_nr[:, :-1]

        x_booms_nr = x_booms_nr.T
        z_booms_nr = z_booms_nr.T

        NAx, NAz = rombus.NeutralAxis(
            A,
            x_booms_nr,
            z_booms_nr)
        assert np.all(np.abs(NAx) <= 0.001), 'Neutral axis is not along the axis for a symmetric shape'
        assert np.all(np.abs(NAz) <= 0.001), 'Neutral axis is not along the axis for a symmetric shape'

    def test_mo_i(self):
        from dse.detailed.Structures.StructureClasses import Beam
        x = np.atleast_2d(np.hstack(
            (np.linspace(-1, 0, 25)[:-1],
             np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25)[:-1])
        ))
        z = np.atleast_2d(np.hstack(
            (np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25)[:-1],
             np.linspace(-1, 0, 25)[:-1]
             )
        ))
        rombus = Beam(
            width=x.T * np.ones(100),
            height=z.T * np.ones(100),
            length=5,
            cross_section=np.vstack((x, z)) * np.ones((100, 1, 1)),
            material='Al/Si',
            fixing_points=np.array([[0], [0]]) * np.ones(100)
        )
        A = np.ones((np.size(x), 100))

        x_booms_nr, z_booms_nr = np.split(np.reshape(rombus.section, (np.size(rombus.y), 2 * np.shape(rombus.x)[0])), 2, 1)
        if np.all(x_booms_nr[:, 0] == x_booms_nr[:, -1]):
            x_booms_nr = x_booms_nr[:, :-1]
            z_booms_nr = x_booms_nr[:, :-1]

        x_booms_nr = x_booms_nr.T
        z_booms_nr = z_booms_nr.T

        Ixx, Izz, Ixz = rombus.MoI(
            A,
            x_booms_nr,
            z_booms_nr
        )

        assert np.all(Ixx - Izz <= 0.001), 'Symmetrical shape does not have equal moments of inertia'
        assert np.all(Ixz <= 0.001), 'Symmetrical shape does not have Ixz = 0'

    def test_stress_calculations(self):
        from dse.detailed.Structures.StructureClasses import Beam, Force
        x = np.hstack(
            (np.linspace(-1, 0, 25)[:-1],
             np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25)[:-1])
        )
        y = np.hstack(
            (np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25)[:-1],
             np.linspace(-1, 0, 25)[:-1]
             )
        )
        rombus = Beam(
            width=x,
            height=y,
            length=5,
            cross_section='constant',
            material=materials['Al/Si']
        )
        A = np.ones((np.size(x), 100))

        dummy_force_1 = Force(
            magnitude=np.array(
                [
                    [0],
                    [1e3],
                    [0]
                ]
            ),
            point_of_application=np.array(
                [
                    [0],
                    [-5],
                    [0]
                ]
            )
        )
        rombus.add_loading(dummy_force_1)

        stress = rombus.StressCalculations(
            boomArea_nr=A
        )

        assert np.all(stress == stress[0]), 'The stress on the booms due to a point load along the length does not ' \
                                            'result in a constant stress'

        rombus.unload()

        dummy_force_2 = Force(
            magnitude=np.array(
                [
                    [0],
                    [0],
                    [1e3]
                ]
            ),
            point_of_application=np.array(
                [
                    [0],
                    [-5],
                    [0]
                ]
            )
        )

        rombus.add_loading(dummy_force_2)

        stress = rombus.StressCalculations(
            boomArea_nr=A
        )

        stress_split = np.split(stress, 4)

        assert np.all((stress_split[0][1:] - np.flip(stress_split[1], 0)[
                                             :-1]) < 0.0001), 'The top booms have non-symmetrical stresses'
        assert np.all((stress_split[2][1:] - np.flip(stress_split[3], 0)[
                                             :-1]) < 0.0001), 'The bottom booms have non-symmetrical stresses'
        assert np.where(stress == np.min(stress))[0] == 24, 'The maximum compressive stress is not at the top'
        assert np.where(stress == np.max(stress))[0] == 72, 'The maximum tensile stress is not at the bottom'
        assert np.max(stress) == np.abs(np.min(stress)), 'The maximum tensile and compressive stresses are not the same'

    def test_internal_stress(self):
        from dse.detailed.Structures.StructureClasses import Beam, Force
        x = np.hstack(
            (np.linspace(-1, 0, 25)[:-1],
             np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25)[:-1])
        )
        y = np.hstack(
            (np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25)[:-1],
             np.linspace(-1, 0, 25)[:-1]
             )
        )
        rombus = Beam(
            width=x,
            height=y,
            length=5,
            cross_section='constant',
            material=materials['Al/Si']
        )

        dummy_force_2 = Force(
            magnitude=np.array(
                [
                    [0],
                    [0],
                    [1e3]
                ]
            ),
            point_of_application=np.array(
                [
                    [0],
                    [-5],
                    [0]
                ]
            )
        )

        rombus.add_loading(dummy_force_2)
        rombus.InternalStress(0, 0, 2)
        assert np.shape(rombus.t) == (np.size(x), 100), 'rombus.t is not an array of the correct shape'
        assert np.all(rombus.t >= 0.001), 'The minimum thickness is smaller than 1mm'
        assert np.all(rombus.sigma <= 250e6/1.5)
