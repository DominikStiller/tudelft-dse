import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv


def xflr_forces(filename, q, b):
    if type(filename) != str:
        raise TypeError(f"Input 1 should be a string, not {type(filename)}")
    if type(q) != float:
        raise TypeError(f"Input 2 should be a float, not {type(q)}")
    if type(b) != float:
        raise TypeError(f"Input 3 should be a float, not {type(b)}")

    # Rewrite the xflr data file into a readable format. ranging from the negative wingtip to closest to zero root
    with open(filename) as csvfile:
        data = csv.reader(csvfile, delimiter=",")
        Data = []
        for row in data:
            Data.append(row)
    while [] in Data:
        Data.remove([])
    indx = Data.index(['Main Wing'])
    Data = Data[indx+1:]

    xflr_data = pd.DataFrame(Data[1:], columns=Data[0], dtype=float)

    # Cut half wing
    cutpoint = int(round(len(xflr_data['Cl']) / 2))
    xflr_data = xflr_data[:cutpoint]


    ## Application of the forces
    # Assume drag acts through the neutral axis
    application = np.zeros((3, len(xflr_data['  y-span'])))   # [m] (3, n) where n is the length of y-span
    application[0] = np.array(xflr_data['XCP'])             # [m] Point of application in the x direction
    application[1] = np.array(xflr_data['  y-span'])          # [m] Point of application in the y direction

    ## Magnitude of the forces
    # Overwrite the forces in y direction if there are any present.
        # [m2] Array of all the surface areas
    S_array = xflr_data['Chord'] * np.hstack((np.array([b-abs(xflr_data['  y-span'][0])]), abs(application[1][1:] - application[1][:-1])))
    mag = np.zeros((3, len(xflr_data['  y-span'])))                          # [N] (3 x n) where n is the length of y_lst
    mag[0] = - (xflr_data['ICd'] + xflr_data['PCd']) * q * S_array       # [N] Forces in x-direction due to drag
    mag[2] = (xflr_data['Cl'] * q * S_array)                               # [N] Forces in z-direction due to lift

    return Force(magnitude=mag, point_of_application=application)


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
            raise TypeError(
                f"Concatenate forces takes as input 1 Force class or a list, not {type(force_2)}"
                f" data types"
            )

    def separate_forces(self):
        mag = np.split(self.F, np.shape(self.F)[1], axis=1)
        point = np.split(self.application, np.shape(self.application)[1], axis=1)
        f_list = list()
        for i in range(len(mag)):
            f_list.append(Force(magnitude=mag[i], point_of_application=point[i]))
        return f_list


class Beam:
    def __init__(self, x, y, z, cross_section):
        # Beam dimensions
        self.x = x  # Width
        self.y = y  # Length
        self.z = z  # Height

        if cross_section == "full":
            self.section = np.ones((len(z), len(x)))
        else:
            self.section = cross_section

        self.coords = np.empty((len(y), len(z), len(x)))
        for i in range(len(y)):
            self.coords[i] = self.section

        # Internal forces
        self.f_loading = np.zeros((len(y), 3, 1))

        # Internal moments
        self.m_loading = np.zeros((len(y), 3, 1))

    def unload(self):
        self.f_loading = np.zeros((len(self.y), 3, 1))
        self.m_loading = np.zeros((len(self.y), 3, 1))

    def add_loading(self,
                    force):
        # Locate where in the beam are the loads located
        for i in range(np.shape(force.F)[1]):
            y_point = (np.abs(self.y - force.application[1, i])).argmin()

            # Add step loads
            self.f_loading[y_point:] += np.reshape(force.F[:, i], (3, 1))

            # Calculate moments
            for j in range(len(self.y[y_point:])):
                distance = np.array([[np.mean(self.x), self.y[y_point+j], np.mean(self.z)]]).T - np.reshape(force.application[:, i], (3, 1))
                self.m_loading[y_point+j] += np.reshape(np.cross(force.F[:, i], distance.T), (3, 1))

    def plot_internal_loading(self):
        fig, (axs1, axs2) = plt.subplots(2, 1)
        f_loading = np.reshape(self.f_loading, (len(self.y), 3))
        m_loading = np.reshape(self.m_loading, (len(self.y), 3))

        axs1.plot(self.y, f_loading[:, 0], label=r'$V_x$')
        axs1.plot(self.y, f_loading[:, 1], label=r'$V_y$')
        axs1.plot(self.y, f_loading[:, 2], label=r'$V_z$')
        axs1.set_xlabel('Span [m]')
        axs1.set_ylabel('Force [N]')
        axs1.legend()

        axs2.plot(self.y, m_loading[:, 0], label=r'$M_x$')
        axs2.plot(self.y, m_loading[:, 1], label=r'$M_y$')
        axs2.plot(self.y, m_loading[:, 2], label=r'$M_z$')
        axs2.set_xlabel('Span [m]')
        axs2.set_ylabel('Moment [Nm]')
        axs2.legend()

        plt.tight_layout()
        plt.show()


if __name__ == '__main__':
    q = 0.5 * 0.01 * 112 ** 2
    aerodynamic_forces = xflr_forces('Test_xflr5_file.csv', q, 16.8)
    x = np.linspace(0, 3.35, 100)
    y = np.linspace(-16.8, 0, 100)
    z = np.linspace(0, 3.35/12, 100)
    wing = Beam(x, y, z, 'full')

    wing.add_loading(aerodynamic_forces)
    wing.plot_internal_loading()
