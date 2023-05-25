import matplotlib.pyplot as plt
import numpy as np
import csv


class Force:
    def __init__(self, magnitude, point_of_application):
        self.F = magnitude
        self.fx = magnitude[0]
        self.fy = magnitude[1]
        self.fz = magnitude[2]
        self.application = point_of_application

    def concatenate_forces(self, force_2):
        if type(force_2) == Force:
            mag = np.hstack((self.F, force_2.F))
            point = np.hstack((self.application, force_2.application))
            return Force(magnitude=mag, point_of_application=point)
        elif type(force_2) == list:
            mag = self.F
            point = self.application
            for force in force_2:
                mag = np.hstack((mag, force.F))
                point = np.hstack((point, force.application))
            return Force(magnitude=mag, point_of_application=point)
        else:
            raise TypeError(f'Concatenate forces takes as input 1 Force class or a list, not {type(force_2)}'
                            f' data types')

    def separate_forces(self):
        mag = np.split(self.F, np.shape(self.F)[1], axis=1)
        point = np.split(self.application, np.shape(self.application)[1], axis=1)
        f_list = list()
        for i in range(len(mag)):
            f_list.append(Force(magnitude=mag[i], point_of_application=point[i]))
        return f_list

    def xflr_forces(self, filename, q):
        # Rewrite the xflr data file into a readable format. ranging from the negative wingtip to closest to zero root
        with open(filename) as csvfile:
            data = csv.reader(csvfile, delimiter=",")
            Data = []
            for row in data:
                Data.append(row)
        y_lst = np.array([float(i[0]) for i in Data[1:]])  # [m] Y-step of the applied loads b
        c_lst = np.array([float(i[1]) for i in Data[1:]])  # [m] Y-step of the applied loads b
        # AoA_lst = [float(i[2]) for i in Data[1:]]  # [degrees] induced angle of angle of attack
        CL_lst = np.array([float(i[3]) for i in Data[1:]])  # [-] Lift coefficient of the wing
        PCd_lst = np.array([float(i[4]) for i in Data[1:]])  # [-] Parasite drag coefficient of the wing
        ICd_lst = np.array([float(i[5]) for i in Data[1:]])  # [-] Induced drag coefficient of the wing
        # CmGeom_lst = [float(i[6]) for i in Data[1:]]  # [-] Moment coefficient
        # CmAirfchord4_lst = [float(i[7]) for i in Data[1:]]  # [-] Moment coefficient at quarter chord
        # XTrtop = [float(i[8]) for i in Data[1:]]  # [m]
        # XTrBot = [float(i[9]) for i in Data[1:]]  # [m]
        XCP = np.array([float(i[10]) for i in Data[1:]])  # [m] x-coordinate of the center of pressure
        # BM = [float(i[11]) for i in Data[1:]]  # [-] Moment coefficient at quarter chord

        ## Application of the forces
        # Assume drag acts through the neutral axis
        self.application = np.zeros((3, len(y_lst)))  # [m] (3, n) where n is the length of y_lst
        self.application[0] = np.array(XCP)  # [m] Point of application in the x direction
        self.application[1] = np.array(y_lst)  # [m] Point of application in the y direction

        ## Magnitude of the forces
        # Overwrite the forces in y direction if there are any present.
        S_array = c_lst*abs(self.application[1][1:]-self.application[1][:-1])  # [m2] Array of all the surface areas
        self.mag = np.zeros((3, len(y_lst))) # [N] (3 x n) where n is the length of y_lst
        self.mag[0][1:] = (ICd_lst+PCd_lst)* q * S_array  # [N] Forces in x-direction due to drag
        self.mag[2] = CL_lst * q * S_array  # [N] Forces in z-direction due to lift






class Beam:
    def __init__(self, x, y, z, cross_section):
        # Beam dimensions
        self.x = x  # Width
        self.y = y  # Length
        self.z = z  # Height

        if cross_section == 'full':
            self.section = np.ones((x, z))
        else:
            self.section = cross_section

        # Internal forces
        self.fx_loading = np.zeros(len(self.x))
        self.fy_loading = np.zeros(len(self.y))
        self.fz_loading = np.zeros(len(self.z))

        # Internal moments
        self.mx_loading = np.zeros(100)
        self.my_loading = np.zeros(100)
        self.mz_loading = np.zeros(100)

    def unload(self):
        self.fx_loading = np.zeros(100)
        self.fy_loading = np.zeros(100)
        self.fz_loading = np.zeros(100)

        self.mx_loading = np.zeros(100)
        self.my_loading = np.zeros(100)
        self.mz_loading = np.zeros(100)

    def add_loading(self, force):
        x_index = (np.abs(self.x - force.application[0])).argmin()
        y_index = (np.abs(self.y - force.application[1])).argmin()
        z_index = (np.abs(self.z - force.application[2])).argmin()

        # Add step loads
        self.fx_loading[:x_index + 1] = force.fx
        self.fz_loading[:y_index + 1] = force.fy
        self.fz_loading[:z_index + 1] = force.fz

        # Calculate moments



