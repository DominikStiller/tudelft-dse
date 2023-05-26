import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import csv
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit


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
    def __init__(self, width, length, height, cross_section):
        # Beam dimensions
        if type(width) == int or type(width) == float:
            self.x = np.reshape(np.linspace(0, width, 100), (100, 1)) * np.ones(100)
        elif type(width) == np.ndarray:
            self.x = width
        else:
            raise TypeError('Width needs to be either a constant float/int or a 100 x n array')

        if type(length) == int or type(length) == float:
            self.y = np.linspace(-length, 0, 100)
        elif type(length) == np.ndarray:
            self.y = length
        else:
            raise TypeError('Length needs to be either a constant float/int or a 100 x n array')

        if type(height) == int or type(height) == float:
            self.z = np.reshape(np.linspace(0, height, 100), (100, 1)) * np.ones(100)
        elif type(height) == np.ndarray:
            self.z = height
        else:
            raise TypeError('Height needs to be either a constant float/int or a 100 x n array')

        if cross_section == "full":
            self.section = np.ones((len(self.y), len(self.z), len(self.x)))
        else:
            if type(cross_section) != np.ndarray:
                raise TypeError('The cross section needs to be either "full" or a 3D array')
            else:
                if np.shape(cross_section) != (len(self.y), len(self.z), len(self.x)):
                    raise ValueError('Cross section needs to have a size consistent with the '
                                     'amount of points indicated for length, height and width')
                self.section = cross_section

        # Internal forces
        self.f_loading = np.zeros((len(self.y), 3, 1))

        # Internal moments
        self.m_loading = np.zeros((len(self.y), 3, 1))

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

    def MoI(self, boom_coordinates, interconnection):
        Airfoil = pd.read_csv('S1223.dat', delimiter='\s+', dtype=float, skiprows=1, names=['x', 'z'])
        if interconnection != 0:
            print('Interconnection between stringers still need to be implemented')
        if boom_coordinates != 0:
            print('The code now only runs for all coordinates of the stringers being the booms')
        else:
            boom_Data = Airfoil

        IndxSplit = np.where(Airfoil['x'] == np.min(Airfoil['x']))

        # Initial gues for the boom areas add up to a percentage fill of the airfoil shape
        # Assume for the first estimate all the boom areas to be equal: Bi = A_airfoil / #booms
        Infill = 0.10  # The percentage infull of the wingbox area

        Airfoil_top = InterpolatedUnivariateSpline(np.flip(Airfoil['x'][:IndxSplit[0][0]+1].to_numpy()),
                                                   np.flip(Airfoil['z'][:IndxSplit[0][0]+1].to_numpy()))
        Airfoil_bot = InterpolatedUnivariateSpline(Airfoil['x'][IndxSplit[0][0]:].to_numpy(),
                                                   Airfoil['z'][IndxSplit[0][0]:].to_numpy())

        A_airfoil = (Airfoil_top.integral(0, 1) - Airfoil_bot.integral(0, 1)) * self.x[-1]  # The Area of the airfoil
        x_booms = np.reshape(boom_Data['x'].to_numpy(), (len(boom_Data['x']), 1)) * self.x[-1]
        z_booms = np.reshape(boom_Data['z'].to_numpy(), (len(boom_Data['z']), 1)) * self.x[-1]

        boomDistance = np.sqrt((x_booms[1:] - x_booms[:-1])**2 + (z_booms[1:] - z_booms[:-1])**2)
        tSkin = (Infill * A_airfoil) / np.sum(boomDistance, 0) * np.ones((np.shape(boomDistance)[0], 1))

        x_booms_nr = x_booms[:-1]
        z_booms_nr = z_booms[:-1]

        Bi_initial = A_airfoil / len(x_booms_nr)
        boomArea_nr = np.ones((np.shape(x_booms_nr)[0], 1)) * Bi_initial * Infill

        diff = np.ones(np.shape(boomArea_nr)[1])
        while np.any(np.abs(diff) > 0.01):
            ## Neutral Axis Calculations
            NAx = (np.sum(boomArea_nr * x_booms_nr, 0)) / (np.sum(boomArea_nr, 0))
            NAz = (np.sum(boomArea_nr * z_booms_nr, 0)) / (np.sum(boomArea_nr, 0))

            Ixx = np.sum(boomArea_nr * (x_booms_nr - NAz) ** 2, 0)
            Izz = np.sum(boomArea_nr * (z_booms_nr - NAx) ** 2, 0)
            Ixz = np.sum(boomArea_nr * (x_booms_nr - NAx) * (z_booms_nr - NAz), 0)

            sigma = ((np.reshape(self.m_loading[:, 0], (np.shape(Izz))) * Izz - np.reshape(self.m_loading[:, 2], (np.shape(Ixz))) * Ixz) * z_booms_nr +
                     (np.reshape(self.m_loading[:, 2], (np.shape(Ixx))) * Ixx - np.reshape(self.m_loading[:, 0], (np.shape(Ixz))) * Ixz) * x_booms_nr) / (Ixx*Izz + Ixz**2)
            sigma = np.vstack((sigma, sigma[0]))
            sigma[np.where(sigma == 0)] = 0.1
            boomAreaCopy = np.copy(boomArea_nr)
            for i in range(np.shape(x_booms_nr)[0]):
                boomAreaCopy[i] = tSkin[i-1] * boomDistance[i-1] / 6 * (2 + sigma[i-1]/sigma[i]) + \
                                  tSkin[i] * boomDistance[i] / 6 * (2 + sigma[i+1]/sigma[i])

            diff = np.mean((boomAreaCopy - boomArea_nr) / boomArea_nr, 0)
            boomArea_nr = np.copy(boomAreaCopy)
            print(np.max(np.abs(diff)), np.where(np.abs(diff) == np.max(np.abs(diff))))
            time.sleep(0.5)



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
    wing = Beam(3.35, 16.8, 3.35/12, 'full')

    wing.add_loading(aerodynamic_forces)
    wing.MoI(0, 0)

    wing.plot_internal_loading()
