import numpy as np


class Fuselage:

    def __init__(self, cabin_diameter, cabin_length, L_D_nose, L_D_tail):
        # Useful info
        self.L_D_nose = L_D_nose
        self.L_D_tail = L_D_tail

        # Useful characteristics
        self.cab_d = cabin_diameter
        self.cab_l = cabin_length

        self.d_main = self.cab_d
        self.r_main = self.d_main / 2

        self.length = self.cab_l + (self.L_D_nose * self.r_main) + (self.L_D_tail * self.r_main)
        r_nose = ((self.L_D_nose * self.r_main) + self.r_main)/2
        self.s = 2 * np.pi * self.r_main * self.cab_l + 0.5 * 4 * np.pi * (
                self.r_main ** 2 + r_nose ** 2)

        self.cd = self.calculate_cd()

        # ---- Legend ---- #
        """
        Height factor describes how much of the total fuselage diameter will not contribute to cabin volume (where the
        astronauts will be at). This "wasted" volume is to be used for wingbox connections.
        
        L_D_nose describes the ratio between the nose's length and its diameter L/D. The same holds for L_D_tail.
        
        cabin_diameter and _length are the useful cabin dimensions (again, where the astronauts at). Those are basically
        the dimensions of the cockpit, the dimensions the astronauts will perceive.
        
        Self d-main and r_main are the effective diameter and radius of the fuselage. So not only useful cockpit 
        dimensions, but the ones seen from the outside and by the flow.
        
        Self length is the total length of the fuselage, from nose tip to tail.
        
        Self S is the external area the flow sees. What we use in Cd*q*S. Not to be confused with S_factor (see 
        explanations inside the function definition of calculate_cd), which Torenbeek lumps together into Cd 
        when he gives us CdS.
        
        Self Cd is the fuselage's Cd.
        """

    #  --------- Methods for calculating characteristics --------- #

    # def effective_fuselage_diameter(self, height_factor):
    #     """ Returns dimensions [...]. Inputs are cabin_diameter, cabin_length, height_factor, all in [m] """
    #     d = self.cab_d * (1 + height_factor)
    #     return d

    def calculate_cd(self):
        """ Calculates cd [-] based on the fuselage characteristics """
        # Calculations in excel - based on Torenbeek

        Cd = 0.0031 * 1 * self.length * 2*self.d_main  # (Torenbeek page 150) # This is Cd * S_factor. Not real S.

        return Cd

    def drag_simulation(self, velocity, air_density):
        """ Calculates fuselage drag [N] based on dynamic pressure. It takes velocity [m/s] and air density [kg/m^3] """
        q = 0.5 * air_density * velocity ** 2
        drag_force = self.cd * self.s * q
        return drag_force


    # --------- Static Methods --------- #

    @staticmethod
    def area_cylinder(d, l):
        """ Returns cylinder's area. Inputs are d = diameter and l = length, both in [m]. Output is in [m^2] """
        r = d / 2
        area = 2 * np.pi * r * (r + l)
        return area

    @staticmethod
    def area_cabin(d, l):
        """ Returns cabin's area. Inputs are d = diameter and l = length, both in [m]. Output is in [m^2] """
        r = d / 2
        area = 2 * np.pi * r * l
        return area


    @staticmethod
    def volume_cylinder(d, l):
        """ Returns cylinder's volume. Inputs are d = diameter and l = length, both in [m]. Output is in [m^3] """
        r = d / 2
        volume = l * np.pi * r ** 2
        return volume

    @staticmethod
    def area_sphere(d):
        """ Returns sphere's area. Input is d = diameter in [m]. Output is in [m^2]"""
        r = d / 2
        area = 4 * np.pi * r ** 2
        return area

    @staticmethod
    def volume_sphere(d):
        """ Returns sphere's volume. Input is d = diameter in [m]. Output is in [m^3]"""
        r = d / 2
        area = (4 / 3) * np.pi * r ** 3
        return area


# Testing

cabin_diameter = 1.66
cabin_length = 2
L_D_nose = 1.5
L_D_tail = 1

fuselage = Fuselage(cabin_diameter, cabin_length, L_D_nose, L_D_tail)

print(f' D:{fuselage.cab_d}        L:{fuselage.length}       S:{round(fuselage.s, 3)}        cd:{round(fuselage.cd, 5)}')

v = 111
rho = 0.015
drag = fuselage.drag_simulation(v, rho)
# drag = round(drag, 3)
print("drag force: ", drag)

"""
#  --------------------------------  TO DO  --------------------------------  #

-> Fix calculations and results based on the new excel sheet
-> Apply changes to the unit tests.

"""
print()

nose_area = Fuselage.area_sphere(2.075)/2
cabin_area = Fuselage.area_cabin(1.66, 2)
rear_area = Fuselage.area_sphere(1.66)/2

total_area = nose_area + cabin_area + rear_area
print(f'Nose area: {nose_area}  \n  cabin area: {cabin_area} \n  rear area: {rear_area}  \n  total area: {total_area}')
