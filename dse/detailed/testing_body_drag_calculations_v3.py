import unittest
from body_drag_calculations_v3 import Fuselage  # Change to V2 later
import numpy as np

"""
Note the name of the tests, e.g. "test_layer1_object_creation".
Layer 1 means how fundamental the functionality that is being testes is, 
in such a way that if a layer 3 test works, all the lower layer tests 
(2 and 1 in this case) are guaranteed to work. So, if a layer 3 fails,
it might no be too hard to fix, but if a layer 1 fails, something really 
fundamental might be broken.
"""

class TestFuselageLogic(unittest.TestCase):

    def test_layer1_object_creation(self):

        cabin_width = 1.5
        cabin_height = 1.8
        fuselage_length = 5.72
        fuselage_area = 25.88


        fus = Fuselage(cabin_width, cabin_height, fuselage_length, fuselage_area)

        self.assertIsInstance(fus, Fuselage, "It appears the class assignment is not working")


    def test_layer3_assign_properties(self):

        cabin_width = 1.5
        cabin_height = 1.8
        fuselage_length = 5.72
        fuselage_area = 25.88

        fus = Fuselage(cabin_width, cabin_height, fuselage_length, fuselage_area)

        assert_equal_message_firsttest_1 = "The initializer values could not be assigned properly"
        assert_equal_message_firsttest_2 = "The fuselage's total length is not being calculated properly"
        assert_equal_message_firsttest_3 = "The fuselage's total surface is not being calculated properly"
        assert_equal_message_firsttest_4 = "The fuselage's Cd is not being calculated properly"


        self.assertEqual(fus.width, cabin_width, assert_equal_message_firsttest_1)
        self.assertEqual(fus.length, fuselage_length, assert_equal_message_firsttest_1)
        self.assertEqual(fus.height, cabin_height, assert_equal_message_firsttest_1)

        self.assertEqual(fus.area, fuselage_area, assert_equal_message_firsttest_3)
        self.assertEqual(fus.cd, fus.calculate_cd(), assert_equal_message_firsttest_4)

        self.assertEqual(fus.width, 1.5, assert_equal_message_firsttest_1)
        self.assertEqual(fus.height, 1.8, assert_equal_message_firsttest_1)
        self.assertEqual(fus.length, 5.72, assert_equal_message_firsttest_1)
        self.assertEqual(round(fus.area, 2), 25.88, assert_equal_message_firsttest_3)
        self.assertEqual(round(fus.cd, 5), 0.02941, assert_equal_message_firsttest_4)




    def test_layer2_check_methods(self):

        cabin_width = 1.5
        cabin_height = 1.8
        fuselage_length = 5.72
        fuselage_area = 25.88

        fus = Fuselage(cabin_width, cabin_height, fuselage_length, fuselage_area)

        self.assertEqual(round(fus.calculate_cd(), 5), 0.02941, "The method calculate_cd is failing to calculate correctly.")


    def test_layer4_check_drag_calculation(self):

        cabin_width = 1.5
        cabin_height = 1.8
        fuselage_length = 5.72
        fuselage_area = 25.88

        fus = Fuselage(cabin_width, cabin_height, fuselage_length, fuselage_area)


        v = 111
        rho = 0.01
        drag = fus.drag_simulation(v, rho)

        self.assertEqual(round(drag, 2), 46.88, "The final drag force calculation does not return the correct value")


    def test_static_methods(self):

        sphere_diameter = 4
        cylinder_diameter = 2
        cylinder_length = 10

        self.assertEqual(Fuselage.area_cylinder(cylinder_diameter, cylinder_length), 22*np.pi, "Test related to static methods failed")
        self.assertEqual(Fuselage.volume_cylinder(cylinder_diameter, cylinder_length), 10*np.pi, "Test related to static methods failed")
        self.assertEqual(Fuselage.area_sphere(sphere_diameter), 16*np.pi, "Test related to static methods failed")
        self.assertEqual(Fuselage.volume_sphere(sphere_diameter), np.pi*32/3, "Test related to static methods failed")


if __name__ == '__main__':
    unittest.main()

