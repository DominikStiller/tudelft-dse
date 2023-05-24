import matplotlib.pyplot as plt
import numpy as np


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


class Structure:
    def __init__(self, structure_type):
        self.type = structure_type


class Beam(Structure):
    def __init__(self, x, y, z, structure_type):
        super().__init__(structure_type)

        # Beam dimensions
        self.x = np.linspace(0, x, 100)
        self.y = np.linspace(0, y, 100)
        self.z = np.linspace(0, z, 100)

        # Internal forces
        self.fx_loading = np.zeros(100)
        self.fy_loading = np.zeros(100)
        self.fz_loading = np.zeros(100)

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



