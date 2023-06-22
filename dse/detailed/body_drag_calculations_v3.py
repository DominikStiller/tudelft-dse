import numpy as np


class Fuselage:
    def __init__(self, cabin_width, cabin_height, fuselage_length, fuselage_area):
        # Useful characteristics
        self.width = cabin_width
        self.height = cabin_height
        self.length = fuselage_length
        self.area = fuselage_area

        self.r_main = (self.width + self.height) / 4
        self.d_main = self.r_main * 2

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

    def calculate_cd(self):
        """Calculates cd [-] based on the fuselage characteristics"""
        # Calculations in excel - based on Torenbeek
        perfect_cylinder_area = 2 * np.pi * self.r_main * (self.r_main + self.length)
        CdS = (
            0.0031 * 1.15 * self.length * (self.height + self.width)
        )  # (Torenbeek page 150) # This is Cd * S_factor. Not real S.
        S_ratio = 0.5 * np.pi * self.length * (self.height + self.width) / 2  # (Torenbeek page 150)

        # The  aerodynamic format of the current shape surely contributes less to drag than a cylinder of equal area.
        # In order to account for this, this ratio is calculated and applied to the hypothetical fully cylindrical
        # fuselage. In summary, it accounts for the better aerodynamics, than just a simple cylinder.

        S_factor = perfect_cylinder_area / S_ratio
        Cd = CdS / S_factor
        return Cd

    def drag_simulation(self, velocity, air_density):
        """Calculates fuselage drag [N] based on dynamic pressure. It takes velocity [m/s] and air density [kg/m^3]"""
        q = 0.5 * air_density * velocity**2
        drag_force = self.cd * self.area * q
        return drag_force

    # --------- Static Methods --------- #

    @staticmethod
    def area_cylinder(d, l):
        """Returns cylinder's area. Inputs are d = diameter and l = length, both in [m]. Output is in [m^2]"""
        r = d / 2
        area = 2 * np.pi * r * (r + l)
        return area

    @staticmethod
    def area_cabin(d, l):
        """Returns cabin's area. Inputs are d = diameter and l = length, both in [m]. Output is in [m^2]"""
        r = d / 2
        area = 2 * np.pi * r * l
        return area

    @staticmethod
    def volume_cylinder(d, l):
        """Returns cylinder's volume. Inputs are d = diameter and l = length, both in [m]. Output is in [m^3]"""
        r = d / 2
        volume = l * np.pi * r**2
        return volume

    @staticmethod
    def area_sphere(d):
        """Returns sphere's area. Input is d = diameter in [m]. Output is in [m^2]"""
        r = d / 2
        area = 4 * np.pi * r**2
        return area

    @staticmethod
    def volume_sphere(d):
        """Returns sphere's volume. Input is d = diameter in [m]. Output is in [m^3]"""
        r = d / 2
        area = (4 / 3) * np.pi * r**3
        return area


if __name__ == "__main__":
    cabin_width = 1.5
    cabin_height = 1.8
    fuselage_length = 6
    fuselage_area = 30.81

    fuselage = Fuselage(cabin_width, cabin_height, fuselage_length, fuselage_area)

    print(
        f" Width:{fuselage.width}        Height:{fuselage.height}       Length:{round(fuselage.length, 3)}    Area:{round(fuselage.area, 5)}     cd:{round(fuselage.cd, 5)}"
    )

    v = 111
    rho = 0.01
    drag = round(fuselage.drag_simulation(v, rho), 2)
    # drag = round(drag, 3)
    print("drag force: ", drag, "[N]")
