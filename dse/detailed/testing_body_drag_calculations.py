import unittest
from body_drag_calculations import Fuselage

class TestFuselageLogic(unittest.TestCase):

    # material.py (Material class)

    def test_layer2_object_creation(self):

        cabin_diameter = 2
        cabin_length = 2
        height_factor = 1/3
        L_D_nose = 1.5
        L_D_tail = 1

        fus = Fuselage(cabin_diameter, cabin_length, height_factor, L_D_nose, L_D_tail)

        self.assertIsInstance(fus, Fuselage, "It appears the class assignment is not working")


    def test_layer1_assign_properties(self):

        cabin_diameter = 2
        cabin_length = 2
        height_factor = 1/3
        L_D_nose = 1.5
        L_D_tail = 1

        fus = Fuselage(cabin_diameter, cabin_length, height_factor, L_D_nose, L_D_tail)

        assert_equal_message_firsttest_1 = "The initializer values could not be assigned properly"
        assert_equal_message_firsttest_2 = "The fuselage's total length is not being calculated properly"
        assert_equal_message_firsttest_3 = "The fuselage's total surface is not being calculated properly"
        assert_equal_message_firsttest_4 = "The fuselage's Cd is not being calculated properly"


        self.assertEqual(fus.cab_d, 2, assert_equal_message_firsttest_1)
        self.assertEqual(fus.d_main, 8/3, assert_equal_message_firsttest_1)
        self.assertEqual(fus.length, 26/3, assert_equal_message_firsttest_2)
        self.assertEqual(round(fus.s, 3), 45.379, assert_equal_message_firsttest_3)
        self.assertEqual(round(fus.cd, 5), 0.03105, assert_equal_message_firsttest_4)


    def test_layer2_check_methods(self):

        cabin_diameter = 2
        cabin_length = 2
        height_factor = 1/3
        L_D_nose = 1.5
        L_D_tail = 1

        fus = Fuselage(cabin_diameter, cabin_length, height_factor, L_D_nose, L_D_tail)

        self.assertEqual(round(fus.calculate_cd(), 5), 0.03105, "The method calculate_cd is failing to calculate correctly.")





        # # Check if the material is indeed an instance of Material.
        #             self.assertIsInstance(material, Material, message_instance)

if __name__ == '__main__':
    unittest.main()
