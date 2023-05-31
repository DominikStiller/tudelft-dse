
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
        self.r_main = self.d_main/2

        self.length = self.cab_l + (self.L_D_nose * self.d_main) + (self.L_D_tail * self.d_main)
        self.s = 2*np.pi*self.r_main * self.cab_l + 0.5*4*np.pi * ( self.r_main**2 + ((self.d_main + self.L_D_nose * self.d_main)/4)**2 )

        self.cd = self.calculate_cd()




    #  --------- Methods for calculating characteristics --------- #

    def effective_cabin_diameter(self, height_factor):
        ''' Returns dimensions [...]. Inputs are cabin_diameter, cabin_length, height_factor, all in [m]'''
        d = self.cab_d * (1 + height_factor)
        return d


    def calculate_cd(self):
        '''Calculates cd based on the fuselage characteristics'''
        # Calculations in excel - based on Torenbeek
        perfect_cylinder_area = 2 * np.pi * self.r_main * (self.r_main + self.length)
        CdS = 0.0031 * 1 * self.length * self.d_main        # (Torenbeek page 150)
        S_ratio = 0.5 * np.pi * self.length * self.d_main   # (Torenbeek page 150)
        S_factor = perfect_cylinder_area / S_ratio
        Cd = CdS / S_factor
        return Cd


    def drag_simulation(self, velocity, air_density):
        ''' Calculates fuselage drag based on dynamic pressure '''
        q = 0.5 * air_density * velocity**2
        drag_force = self.cd * self.s * q
        return drag_force


    # --------- Static Methods --------- #

    @staticmethod
    def area_cylinder(d,l):
        ''' Returns cylinder area. Inputs are d = diameter and l = length, both in [m]'''
        r = d/2
        area = 2*np.pi*r * (r + l)
        return area


    @staticmethod
    def volume_cylinder(d,l):
        ''' Returns cylinder volume. Inputs are d = diameter and l = length, both in [m]'''

        r = d/2
        volume = l*np.pi*r**2
        return volume


    @staticmethod
    def area_sphere(d):
        ''' Returns sphere area. Input is d = diameter in [m]'''

        r = d/2
        area = 4*np.pi*r**2
        return area


    @staticmethod
    def volume_sphere(d):
        ''' Returns sphere area. Input is d = diameter in [m]'''

        r = d/2
        area = (4/3)*np.pi*r**3
        return area


# Testing

fuselage = Fuselage(2, 2, 1/3, 1.5, 1)

print(f'{fuselage.cab_d}         {fuselage.d_main}         {fuselage.length}        {fuselage.cd}')




"""
#  --------------------------------  TO DO  --------------------------------  #

-> Document assumptions

-> Write proper class descriptions

-> Format it to be a pretty code


"""



