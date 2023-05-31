import unittest
from body_drag_calculations import Fuselage

class TestFuselageLogic(unittest.TestCase):

    # material.py (Material class)

    def test_layer1_object_creation(self):

        cabin_diameter = 8/3
        cabin_length = 2
        L_D_nose = 1.5
        L_D_tail = 1

        fus = Fuselage(cabin_diameter, cabin_length, L_D_nose, L_D_tail)

        self.assertIsInstance(fus, Fuselage, "It appears the class assignment is not working")


    def test_layer3_assign_properties(self):

        cabin_diameter = 8/3
        cabin_length = 2
        L_D_nose = 1.5
        L_D_tail = 1

        fus = Fuselage(cabin_diameter, cabin_length, L_D_nose, L_D_tail)

        assert_equal_message_firsttest_1 = "The initializer values could not be assigned properly"
        assert_equal_message_firsttest_2 = "The fuselage's total length is not being calculated properly"
        assert_equal_message_firsttest_3 = "The fuselage's total surface is not being calculated properly"
        assert_equal_message_firsttest_4 = "The fuselage's Cd is not being calculated properly"


        self.assertEqual(fus.cab_d, 8/3, assert_equal_message_firsttest_1)
        self.assertEqual(fus.length, 26/3, assert_equal_message_firsttest_2)
        self.assertEqual(round(fus.s, 3), 45.379, assert_equal_message_firsttest_3)
        self.assertEqual(round(fus.cd, 5), 0.03105, assert_equal_message_firsttest_4)


    def test_layer2_check_methods(self):

        cabin_diameter = 8/3
        cabin_length = 2
        L_D_nose = 1.5
        L_D_tail = 1

        fus = Fuselage(cabin_diameter, cabin_length, L_D_nose, L_D_tail)

        self.assertEqual(round(fus.calculate_cd(), 5), 0.03105, "The method calculate_cd is failing to calculate correctly.")


    def test_layer4_check_drag_calculation(self):

        cabin_diameter = 8/3
        cabin_length = 2
        L_D_nose = 1.5
        L_D_tail = 1

        fus = Fuselage(cabin_diameter, cabin_length, L_D_nose, L_D_tail)

        v = 111
        rho = 0.015
        drag = fus.drag_simulation(v, rho)

        self.assertEqual(round(drag, 3), 130.185, "The final drag force calculation does not return the correct value")


if __name__ == '__main__':
    unittest.main()


"""
#  --------------------------------  TO DO  --------------------------------  #

-> Implement tests for static methods


"""