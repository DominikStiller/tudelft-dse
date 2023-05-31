from unittest import TestCase
import numpy as np
import os


class Test(TestCase):
    def test_xflr_forces(self):
        from dse.detailed.Structures.loading import xflr_forces

        if os.getcwd().split("\\")[-1] == "structures":
            os.chdir("..\\..\\..\\dse\\detailed\\Structures")
        elif os.getcwd().split("\\")[-1] != "Structures":
            os.chdir("..\\dse\\detailed\\Structures")

        force = xflr_forces(filename="Test_xflr5_file.csv", q=0.5 * 0.01 * 112**2, b=16.8)

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
        from dse.detailed.Structures.loading import xflr_forces

        try:
            Forces = xflr_forces(5, q=0.5 * 0.01 * 112**2, b=16.8)
            a = False
        except TypeError:
            a = True
        self.assertTrue(a)

    def test_xflr_forces_q_type(self):
        from dse.detailed.Structures.loading import xflr_forces

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
        from dse.detailed.Structures.loading import xflr_forces

        if os.getcwd().split("\\")[-1] == "structures":
            os.chdir("..\\..\\..\\dse\\detailed\\Structures")
        elif os.getcwd().split("\\")[-1] != "Structures":
            os.chdir("..\\dse\\detailed\\Structures")

        try:
            Forces = xflr_forces(filename="Test_xflr5_file.csv", q=0.5 * 0.01 * 112**2, b="16.8")
            a = False
        except TypeError:
            a = True
        self.assertTrue(a)


class TestForce(TestCase):
    def test_concatenate_two_forces(self):
        from dse.detailed.Structures.loading import Force

        dummy = Force(
            magnitude=np.array([[3, 3, 3]]).T, point_of_application=np.array([[1, 0, 0]]).T
        )
        dummy_2 = Force(
            magnitude=np.array([[-2, 0, 1]]).T, point_of_application=np.array([[2, 0, 0]]).T
        )
        dummy = dummy.concatenate_forces(dummy_2)

        assert np.shape(dummy.F) == (3, 2)

    def test_concatenate_functions_list(self):
        from dse.detailed.Structures.loading import Force

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
        from dse.detailed.Structures.loading import Force

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
        from dse.detailed.Structures.loading import Beam, Force

        x = np.linspace(0, 0.1, 10)
        y = np.linspace(-2, 0, 100)
        z = np.linspace(0, 0.05, 5)
        wing = Beam(x, y, z, "full")

        F = 1
        loc = -1.5

        point_load = Force(
            magnitude=np.array([[0, 0, F]]).T, point_of_application=np.array([[0, loc, 0]]).T
        )

        wing.add_loading(point_load)

        y_index = (np.abs(y - loc)).argmin()
        if np.all(wing.f_loading[y_index:] == np.array([[0, 0, F]]).T):
            if (np.abs((wing.m_loading[y_index] - wing.m_loading[-1]) / loc)[0] - F) < 0.01:
                assert True
            else:
                assert (np.abs((wing.m_loading[y_index] - wing.m_loading[-1]) / loc)[0] - F) < 0.01
        else:
            assert np.all(wing.f_loading[y_index:] == np.array([[0, 1, 0]]).T)

    def test_add_loading_distributed_load(self):
        from dse.detailed.Structures.loading import Beam, Force

        x = np.linspace(0, 0.1, 10)
        y = np.linspace(-2, 0, 100)
        z = np.linspace(0, 0.05, 5)
        wing = Beam(x, y, z, "full")

        F = np.array([[0, 0, 0], [0, 0, 0], [1, 1, 1]])

        app = np.array([[0.05, 0.05, 0.05], [-0.1, -1, -1.5], [0.025, 0.025, 0.025]])

        force = Force(magnitude=F, point_of_application=app)
        wing.add_loading(force)

        y_index_1 = (np.abs(y - app[:, 0][1])).argmin()
        y_index_2 = (np.abs(y - app[:, 1][1])).argmin()
        y_index_3 = (np.abs(y - app[:, 2][1])).argmin()

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
