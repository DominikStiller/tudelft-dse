from unittest import TestCase
import numpy as np


class Test(TestCase):
    def test_xflr_forces(self):
        from dse.detailed.Structures.loading import xflr_forces

        force = xflr_forces(filename='..\\dse\\detailed\\Structures\\Test_xflr5_file.csv',
                            q=0.5*0.01*112**2,
                            b=16.8)

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
            Forces = xflr_forces(5, q=0.5 * 0.01 * 112 ** 2, b = 16.8)
            a = False
        except TypeError:
            a = True
        self.assertTrue(a)

    def test_xflr_forces_q_type(self):
        from dse.detailed.Structures.loading import xflr_forces
        try:
            Forces = xflr_forces('..\\dse\\detailed\\Structures\\Test_xflr5_file.csv',
                                 q='0.5 * 0.01 * 112 ** 2',
                                 b=16.8)
            a = False
        except TypeError:
            a = True
        self.assertTrue(a)
    def test_xflr_forces_b_type(self):
        from dse.detailed.Structures.loading import xflr_forces

        try:
            Forces = xflr_forces(filename='..\\dse\\detailed\\Structures\\Test_xflr5_file.csv',
                                 q=0.5 * 0.01 * 112 ** 2,
                                 b='16.8')
            a = False
        except TypeError:
            a = True
        self.assertTrue(a)





class TestForce(TestCase):
    def test_concatenate_two_forces(self):
        from dse.detailed.Structures.loading import Force
        dummy = Force(magnitude=np.array([[3, 3, 3]]).T, point_of_application=np.array([[1, 0, 0]]).T)
        dummy_2 = Force(magnitude=np.array([[-2, 0, 1]]).T, point_of_application=np.array([[2, 0, 0]]).T)
        dummy = dummy.concatenate_forces(dummy_2)

        assert np.shape(dummy.F) == (3, 2)

    def test_concatenate_functions_list(self):
        from dse.detailed.Structures.loading import Force
        dummy = Force(magnitude=np.array([[3, 3, 3]]).T,
                      point_of_application=np.array([[1, 0, 0]]).T)
        dummy_2 = Force(magnitude=np.array([[-2, 0, 1]]).T,
                        point_of_application=np.array([[2, 0, 0]]).T)
        dummy_3 = Force(magnitude=np.array([[1, 0, 0]]).T,
                        point_of_application=np.array([[3, 0, 0]]).T)

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
    def test_add_loading(self):
        self.fail()
