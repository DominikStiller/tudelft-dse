import numpy as np


class Fuselage:

    def __init__(self, cabin_diameter, cabin_length, height_factor, L_D_nose, L_D_tail):
        # Useful info
        self.L_D_nose = L_D_nose
        self.L_D_tail = L_D_tail
        self.height_factor = height_factor

        # Useful characteristics
        self.cab_d = cabin_diameter
        self.cab_l = cabin_length

        self.d_main = self.effective_cabin_diameter(self.height_factor)
        self.r_main = self.d_main / 2

        self.length = self.cab_l + (self.L_D_nose * self.d_main) + (self.L_D_tail * self.d_main)
        self.s = 2 * np.pi * self.r_main * self.cab_l + 0.5 * 4 * np.pi * (
                    self.r_main ** 2 + ((self.d_main + self.L_D_nose * self.d_main) / 4) ** 2)

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

    def effective_cabin_diameter(self, height_factor):
        """ Returns dimensions [...]. Inputs are cabin_diameter, cabin_length, height_factor, all in [m] """
        d = self.cab_d * (1 + height_factor)
        return d

    def calculate_cd(self):
        """ Calculates cd [-] based on the fuselage characteristics """
        # Calculations in excel - based on Torenbeek
        perfect_cylinder_area = 2 * np.pi * self.r_main * (self.r_main + self.length)
        CdS = 0.0031 * 1 * self.length * self.d_main  # (Torenbeek page 150) # This is Cd * S_factor. Not real S.
        S_ratio = 0.5 * np.pi * self.length * self.d_main  # (Torenbeek page 150)

        # The  aerodynamic format of the current shape surely contributes less to drag than a cylinder of equal area.
        # In order to account for this, this ratio is calculated and applied to the hypothetical fully cylindrical
        # fuselage. In summary, it accounts for the better aerodynamics, than just a simple cylinder.

        S_factor = perfect_cylinder_area / S_ratio
        Cd = CdS / S_factor
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

fuselage = Fuselage(2, 2, 1 / 3, 1.5, 1)

print(f'{fuselage.cab_d}         {fuselage.d_main}         {fuselage.length}        {fuselage.cd}')

"""
#  --------------------------------  TO DO  --------------------------------  #

V -> Document assumptions

V -> Write proper class descriptions

V -> Format it to be a pretty code


"""
