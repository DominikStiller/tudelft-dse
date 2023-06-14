from dse.detailed.Structures.material_properties import materials, Material
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib as mpl
import vibration_toolbox as vtb
import numpy as np
import pandas as pd
import csv
import scipy as sp
from tqdm import tqdm
from colorama import Fore


def xflr_forces(filename, q, b, adrian=None):
    if type(filename) != str:
        raise TypeError(f"Input 1 should be a string, not {type(filename)}")
    if type(q) != float and type(q) != int:
        raise TypeError(f"Input 2 should be a float or int, not {type(q)}")
    if type(b) != float and type(b) != int:
        raise TypeError(f"Input 3 should be a float or int, not {type(b)}")

    # Rewrite the xflr data file into a readable format. ranging from the negative wingtip to closest to zero root
    if type(adrian) != int:
        with open(filename) as csvfile:
            data = csv.reader(csvfile, delimiter=",")
            Data = []
            for row in data:
                Data.append(row)
        while [] in Data:
            Data.remove([])
        indx = Data.index(["Main Wing"])
        Data = Data[indx + 1 :]

        xflr_data = pd.DataFrame(Data[1:], columns=Data[0], dtype=float)

        # Cut half wing
        cutpoint = int(round(len(xflr_data["Cl"]) / 2))
        xflr_data = xflr_data[:cutpoint]

        ## Application of the forces
        # Assume drag acts through the neutral axis
        application = np.zeros(
            (3, len(xflr_data["  y-span"]))
        )  # [m] (3, n) where n is the length of y-span
        application[0] = np.array(xflr_data["XCP"])  # [m] Point of application in the x direction
        application[1] = np.array(
            xflr_data["  y-span"]
        )  # [m] Point of application in the y direction

        ## Magnitude of the forces
        # Overwrite the forces in y direction if there are any present.
        # [m2] Array of all the surface areas
        S_array = xflr_data["Chord"] * np.hstack(
            (
                np.array([b - abs(xflr_data["  y-span"][0])]),
                abs(application[1][1:] - application[1][:-1]),
            )
        )

        margin = 500
        cd = xflr_data["ICd"] + xflr_data["PCd"]
        lift_drag_ratio = xflr_data["Cl"] / cd
        mag = np.zeros(
            (3, len(xflr_data["  y-span"]))
        )  # [N] (3 x n) where n is the length of y_lst
        mag[0] = -cd * q * S_array + margin / lift_drag_ratio * cd / np.sum(
            cd
        )  # [N] Forces in x-direction due to drag
        mag[2] = xflr_data["Cl"] * q * S_array + margin * xflr_data["Cl"] / np.sum(
            xflr_data["Cl"]
        )  # [N] Forces in z-direction due to lift

        return Force(magnitude=mag, point_of_application=application)

    else:
        if adrian == -1:
            df = pd.read_csv(filename)
            cl_array = np.fromstring(df["Cl_dist"][0][1:-1], sep=" ")
            cd_array = np.fromstring(df["Cd_dist"][0][1:-1], sep=" ")
            sp_array = np.fromstring(df["y_span"][0][1:-1], sep=" ")
        else:
            df = pd.read_csv(filename)
            cl_array = np.fromstring(df["Cl_dist"].iloc[adrian][1:-1], sep=" ")
            cd_array = np.fromstring(df["Cd_dist"].iloc[adrian][1:-1], sep=" ")
            sp_array = np.fromstring(df["y_span"].iloc[adrian][1:-1], sep=" ")

        return cl_array, cd_array, sp_array


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
    def __init__(self, width, length, height, cross_section, material, fixing_points):
        # Beam dimensions
        if type(length) == int or type(length) == float:
            self.y = np.linspace(-length, 0, 100)
        elif type(length) == np.ndarray:
            self.y = length
        else:
            raise TypeError("Length needs to be either a constant float/int or a 1D array")

        if type(width) == int or type(width) == float:
            if cross_section == "square":
                self.x = np.atleast_2d(
                    np.hstack(
                        (
                            np.zeros(20),
                            np.linspace(0, width, 21)[:-1],
                            np.ones(20) * width,
                            np.linspace(width, 0, 21)[:-1],
                        )
                    )
                ).T * np.ones(np.size(self.y))
            else:
                self.x = np.reshape(np.linspace(0, width, 100), (100, 1)) * np.ones(np.size(self.y))
        elif type(width) == np.ndarray:
            if len(np.shape(width)) == 2 and np.shape(width)[1] == np.size(self.y):
                self.x = width
            else:
                raise ValueError("Width array must be n x len(y)")
        else:
            raise TypeError("Width needs to be either a constant float/int or a n x len(y) array")

        if type(height) == int or type(height) == float:
            if cross_section == "square":
                self.z = np.atleast_2d(
                    np.hstack(
                        (
                            np.linspace(0, height, 21)[:-1],
                            np.ones(20) * height,
                            np.linspace(height, 0, 21)[:-1],
                            np.zeros(20),
                        )
                    )
                ).T * np.ones(np.size(self.y))
            else:
                self.z = np.reshape(np.linspace(0, height, 100), (100, 1)) * np.ones(
                    np.size(self.y)
                )
        elif type(height) == np.ndarray:
            if len(np.shape(height)) == 2 and np.shape(height)[1] == np.size(self.y):
                self.z = height
            else:
                raise ValueError("Height array must be n x len(y)")
        else:
            raise TypeError("Height needs to be either a constant float/int or a n x len(y) array")

        if np.shape(self.x) != np.shape(self.z):
            raise ValueError("Width and height must have matching shapes")

        if type(cross_section) == str:
            if cross_section == "constant":
                self.section = np.ones((len(self.y), 2, len(self.x))) * np.vstack((self.x, self.z))
            elif cross_section == "square":
                self.section = np.vstack((self.x[:, 0], self.z[:, 0])) * np.ones(
                    (np.size(self.y), 2, np.shape(self.x)[0])
                )
            else:
                raise ValueError('The cross section needs to be either "constant" or a 3D array')
        else:
            if type(cross_section) != np.ndarray:
                raise TypeError('The cross section needs to be either "constant" or a 3D array')
            else:
                if np.shape(cross_section) != (len(self.y), 2, len(self.x)):
                    raise ValueError(
                        "Cross section needs to have a size (len(self.y), 2, len(self.x))"
                    )
                self.section = cross_section

        # Internal forces
        self.f_loading = np.zeros((len(self.y), 3, 1))

        # Internal moments
        self.m_loading = np.zeros((len(self.y), 3, 1))
        mat = dict()
        self.material_types = dict()
        for i, k in enumerate(materials):
            mat[k] = i
            self.material_types[i] = materials[k]

        if type(material) == str:
            m = mat[material]
            self.mat = m * np.ones(
                (np.shape(self.section)[2], np.shape(self.section)[0]), dtype=int
            )
        elif type(material) == np.ndarray:
            if np.shape(material) == (np.shape(self.section)[2], np.shape(self.section)[0]):
                self.mat = material
            else:
                raise TypeError(
                    "If material is an array, it must be the same shape as the cross-section"
                )
        else:
            raise TypeError("Material needs to be a single Material class or an array")

        if np.shape(fixing_points) == (2, 1) and cross_section == "constant":
            self.fix = fixing_points * np.ones((np.size(self.y), 2, 1))
        elif np.shape(fixing_points) == (2, np.size(self.y)):
            self.fix = fixing_points
        else:
            raise ValueError(
                "Fixing points needs to be a 2x1 array with constant cross section or 2xn array"
            )

        # Parameters to be calculated later
        self.t = None
        self.sigma = None
        self.Bi = None
        self.m = None
        self.Ix = None
        self.Iy = None
        self.Iz = None
        self.xcg = None
        self.ycg = None
        self.zcg = None

    def unload(self):
        self.f_loading = np.zeros((len(self.y), 3, 1))
        self.m_loading = np.zeros((len(self.y), 3, 1))

    def add_loading(self, force):
        # Locate where in the beam are the loads located
        for i in range(np.shape(force.F)[1]):
            y_point = (np.abs(self.y - force.application[1, i])).argmin()

            # Add step loads
            self.f_loading[y_point:] += np.reshape(force.F[:, i], (3, 1))

            # Calculate moments
            for j in range(len(self.y[y_point:])):
                distance = np.array(
                    [[self.fix[0][y_point + j], self.y[y_point + j], self.fix[1][y_point + j]]]
                ).T - np.reshape(force.application[:, i], (3, 1))
                self.m_loading[y_point + j] += np.reshape(
                    np.cross(force.F[:, i], distance.T), (3, 1)
                )

    def add_moment(self, arr, y):
        for i in range(np.shape(arr)[1]):
            y_point = (np.abs(self.y - y)).argmin()

            # Add step loads
            self.m_loading[y_point:] += np.reshape(arr, (3, 1))

    def AirfoilArea(self):
        # Locating the leading edge.
        if len(np.shape(self.x)) == 1:
            IndxSplit = np.where(self.x == np.min(self.x))
            # Calculating the area for the top skin and the bottom skin
            Airfoil_top = InterpolatedUnivariateSpline(
                np.flip(self.x[: IndxSplit[0][0] + 1]), np.flip(self.z[: IndxSplit[0][0] + 1])
            )
            Airfoil_bot = InterpolatedUnivariateSpline(
                self.x[IndxSplit[0][0] :], self.z[IndxSplit[0][0] :]
            )
            A_airfoil = (Airfoil_top.integral(0, 1) - Airfoil_bot.integral(0, 1)) * max(
                self.x
            )  # The Area of the airfoil
        elif len(np.shape(self.x)) == 2:
            A_airfoil = np.zeros(np.shape(self.x)[1])
            for i in range(np.shape(self.x)[1]):
                IndxSplit = np.where(self.x[:, i] == np.min(self.x[:, i]))
                # Calculating the area for the top skin and the bottom skin
                Airfoil_top = InterpolatedUnivariateSpline(
                    np.flip(self.x[:, i][: IndxSplit[0][0] + 1]),
                    np.flip(self.z[:, i][: IndxSplit[0][0] + 1]),
                )
                Airfoil_bot = InterpolatedUnivariateSpline(
                    self.x[:, i][IndxSplit[0][0] :], self.z[:, i][IndxSplit[0][0] :]
                )
                A_airfoil[i] = (Airfoil_top.integral(0, 1) - Airfoil_bot.integral(0, 1)) * max(
                    self.x[:, i]
                )  # The Area of the airfoil
        else:
            raise TypeError("self.x is not a 1D or 2D array")

        return A_airfoil

    def InitialBoom(self, i):
        Infill = 0.10
        ## Reshapes the x and z coordinates to the correct shapes, one with repeating first node, and one without (_nr)
        x_booms = np.reshape(self.x, (len(self.x), 1))
        z_booms = np.reshape(self.z, (len(self.z), 1))
        x_booms_nr = x_booms[:-1]
        z_booms_nr = z_booms[:-1]
        if i == 0:
            A_airfoil = self.AirfoilArea()
        else:
            A_airfoil = 1
        Bi_initial = A_airfoil / len(x_booms_nr) * Infill
        return Bi_initial

    def DesignConstraints(self):
        n = 1.5  # [-] Safety factor
        sigma_yield = np.zeros(np.shape(self.mat))
        for i in range(np.shape(self.mat)[0]):
            for j in range(np.shape(self.mat)[1]):
                sigma_yield[i, j] = self.material_types[
                    self.mat[i, j]
                ].compressive  # [MPa] The yield strength
        tSkin_min = 0.001  # [m] The minimal allowable thickness
        return n, sigma_yield, tSkin_min

    def NeutralAxis(self, boomArea_nr, x_booms_nr, z_booms_nr):
        # Neutral Axis Calculations
        NAx = (np.sum(boomArea_nr * x_booms_nr, 0)) / (np.sum(boomArea_nr, 0))  # [m]
        NAz = (np.sum(boomArea_nr * z_booms_nr, 0)) / (np.sum(boomArea_nr, 0))  # [m]
        return NAx, NAz

    def MoI(self, boomArea_nr, x_booms_nr, z_booms_nr):
        # Moment of Inertia calculations
        Ixx = np.sum(boomArea_nr * (z_booms_nr - self.fix[1]) ** 2, 0)
        Izz = np.sum(boomArea_nr * (x_booms_nr - self.fix[0]) ** 2, 0)
        Ixz = np.sum(boomArea_nr * (x_booms_nr - self.fix[0]) * (z_booms_nr - self.fix[1]), 0)
        return Ixx, Izz, Ixz

    def StressCalculations(self, boomArea_nr):
        Mx = np.reshape(self.m_loading[:, 0], np.size(self.y))
        Mz = np.reshape(self.m_loading[:, 2], np.size(self.y))
        Fy = np.reshape(self.f_loading[:, 1], np.size(self.y))

        x_booms, z_booms = np.split(
            np.reshape(self.section, (np.size(self.y), 2 * np.shape(self.x)[0])), 2, 1
        )
        if np.all(x_booms[:, 0] == x_booms[:, -1]):
            x_booms = x_booms[:, :-1]
            z_booms = z_booms[:, :-1]

        x_booms = x_booms.T
        z_booms = z_booms.T

        # Call the neutral & moments of inertia
        NAx, NAz = self.fix
        Ixx, Izz, Ixz = self.MoI(boomArea_nr, x_booms, z_booms)
        # Stress calculations
        sigma_nr = (
            (Mx * Izz - Mz * Ixz) * (z_booms - NAz) + (Mz * Ixx - Mx * Ixz) * (x_booms - NAx)
        ) / (Ixx * Izz - Ixz**2) - (Fy / np.sum(boomArea_nr, 0)) * np.ones(np.shape(boomArea_nr))
        return sigma_nr

    def TorsionStress(self, boomArea_nr, A_Airfoil):
        Vx = np.reshape(self.f_loading[:, 0], np.size(self.y))
        Vz = np.reshape(self.f_loading[:, 2], np.size(self.y))
        My = np.reshape(self.m_loading[:, 1], np.size(self.y))
        x0 = self.fix[0]
        z0 = self.fix[1]
        if A_Airfoil == 0:
            A = self.AirfoilArea()
        else:
            A = A_Airfoil
        x_booms, z_booms = np.split(
            np.reshape(self.section, (np.size(self.y), 2 * np.shape(self.x)[0])), 2, 1
        )
        if np.all(x_booms[:, 0] == x_booms[:, -1]):
            x_booms_nr = x_booms[:, :-1]
            z_booms_nr = z_booms[:, :-1]
        else:
            x_booms_nr = np.copy(x_booms)
            z_booms_nr = np.copy(z_booms)
            x_booms = np.hstack(
                (x_booms_nr, np.reshape(x_booms_nr[:, 0], (np.shape(x_booms_nr)[0], 1)))
            )
            z_booms = np.hstack(
                (z_booms_nr, np.reshape(z_booms_nr[:, 0], (np.shape(z_booms_nr)[0], 1)))
            )

        x_booms_nr = x_booms_nr.T
        z_booms_nr = z_booms_nr.T
        x_booms = x_booms.T
        z_booms = z_booms.T

        # Call the neutral & moments of inertia
        NAx, NAz = self.NeutralAxis(boomArea_nr, x_booms_nr, z_booms_nr)
        Ixx, Izz, Ixz = self.MoI(boomArea_nr, x_booms_nr, z_booms_nr)

        C = -(Vz * Izz - Vx * Ixz) / (Ixx * Izz - Ixz**2)
        D = -(Vx * Ixx - Vz * Ixz) / (Ixx * Izz - Ixz**2)

        qb = np.zeros(np.shape(boomArea_nr))

        for i in range(1, len(qb[:, 0])):
            qb[i] = (
                qb[i - 1]
                + C * boomArea_nr[i] * (z_booms_nr[i] - NAz)
                + D * boomArea_nr[i] * (x_booms_nr[i] - NAx)
            )
        l = np.sqrt((x_booms[1:] - x_booms[:-1]) ** 2 + (z_booms[1:] - z_booms[:-1]) ** 2)
        d = []

        for i in range(len(x_booms_nr)):
            if x_booms[i, 0] == x_booms[i + 1, 0]:
                d.append((x_booms[i] - x0))
            else:
                a = (z_booms[i + 1] - z_booms[i]) / (x_booms[i] - x_booms[i + 1])
                b = 1
                c = (x_booms[i + 1] * z_booms[i] - x_booms[i] * z_booms[i + 1]) / (
                    x_booms[i] - x_booms[i + 1]
                )
                d.append(abs(a * x0 + b * z0 + c) / np.sqrt(a**2 + b**2))
        d = np.array(d)
        qs0 = (-My - np.sum(qb * l * d, 0)) / (2 * A)
        q_total = qb + np.ones(np.shape(qb)) * qs0
        return q_total

    def BoomArea(self, boomAreaCopy, tSkin, boomDistance, sigma):
        # Update the boom areas for the given stress.
        for i in range(np.shape(boomAreaCopy)[0]):
            boomAreaCopy[i] = tSkin[i - 1] * boomDistance[i - 1] / 6 * (
                2 + sigma[i - 1] / sigma[i]
            ) + tSkin[i] * boomDistance[i] / 6 * (2 + sigma[i + 1] / sigma[i])
        return boomAreaCopy

    def InternalStress(self, boom_coordinates, interconnection, i, title=None):
        if interconnection != 0:  # define the interconnection of all of the boom areas
            print("Interconnection between stringers still need to be implemented")

        if boom_coordinates != 0:  # Defines locations of the booms
            print("The code now only runs for all coordinates of the stringers being the booms")

        x_booms, z_booms = np.split(
            np.reshape(self.section, (np.size(self.y), 2 * np.shape(self.x)[0])), 2, 1
        )
        if np.all(x_booms[:, 0] == x_booms[:, -1]):
            x_booms_nr = x_booms[:, :-1]
            z_booms_nr = z_booms[:, :-1]
        else:
            x_booms_nr = np.copy(x_booms)
            z_booms_nr = np.copy(z_booms)
            x_booms = np.hstack((x_booms, np.reshape(x_booms[:, 0], (np.size(self.y), 1))))
            z_booms = np.hstack((z_booms, np.reshape(z_booms[:, 0], (np.size(self.y), 1))))

        x_booms = x_booms.T
        x_booms_nr = x_booms_nr.T
        z_booms = z_booms.T
        z_booms_nr = z_booms_nr.T

        ## Calculate the distance between each boom, there needs to be an update for the vertical connections
        boomDistance = np.sqrt(
            (x_booms[1:] - x_booms[:-1]) ** 2 + (z_booms[1:] - z_booms[:-1]) ** 2
        )

        ## The design constraints
        n, sigma_ult, tSkin_min = self.DesignConstraints()

        ## Initial guess for the boom area and skin thickness and stress
        tSkin = np.ones(np.shape(boomDistance)) * tSkin_min
        Bi_initial = (
            np.sum(boomDistance, 0)
            * tSkin_min
            / np.shape(boomDistance)[0]
            * np.ones(np.shape(boomDistance))
        )
        boomArea_nr = np.ones((np.shape(x_booms_nr)[0], np.size(self.y))) * Bi_initial
        sigma = n * sigma_ult

        br = False
        while np.any(np.abs(sigma - sigma_ult / n) / (sigma_ult / n) > 0.1):
            diff = np.ones(np.shape(boomArea_nr)[1])
            while np.any(np.abs(diff) > 0.01):
                # Stress Calculations
                sigma_nr = self.StressCalculations(boomArea_nr)

                # Stress with repeating column for the first node for easy calculations for the boom areas
                sigma = np.vstack((sigma_nr, sigma_nr[0]))
                sigma[np.where(np.abs(sigma) <= 0.01)] = 0.1
                q_total = self.TorsionStress(boomArea_nr, i)
                tau = q_total / tSkin
                if np.any(np.abs(sigma) > sigma_ult / n):
                    # Boom area calculations
                    boomAreaCopy = np.copy(boomArea_nr)
                    boomAreaCopy = self.BoomArea(boomAreaCopy, tSkin, boomDistance, sigma)

                    # Calculating the difference between the new and old boom areas.
                    diff = np.mean((boomAreaCopy - boomArea_nr) / boomArea_nr, 0)
                    boomArea_nr = np.copy(boomAreaCopy)
                else:
                    br = True
                    break

            if br:
                break
            indx = list()
            for k in range(np.shape(sigma)[1]):
                if np.any(np.max(abs(sigma[:, k])) >= sigma_ult):
                    indx.append(k)
            if len(indx) > 0:
                for col in indx:
                    tSkin[:, col] += 0.001
            else:
                break
            if np.any(np.isnan(sigma)):
                raise ValueError("There was a NaN in the stress array")

        self.t = tSkin
        self.sigma = sigma_nr
        self.tau = tau
        self.Bi = boomArea_nr
        # plt.axhline(sigma_ult)
        # plt.axhline(-sigma_ult)
        plt.plot(self.y, sigma.transpose())
        if title is not None:
            plt.title(title)
        plt.show()

    def calculate_mass(self):
        dy = self.y[1:] - self.y[:-1]
        volume = dy * self.Bi[:, :-1]
        rho = np.zeros(np.shape(volume))
        for i in range(np.shape(self.mat)[0] - 1):
            for j in range(np.shape(self.mat)[1] - 1):
                rho[i, j] = self.material_types[self.mat[i, j]].rho
        self.m = np.sum(volume * rho)

    def rho(self):
        rho = np.zeros(np.shape((self.Bi[:, :-1] * self.y[1:] - self.y[:-1])))
        if self.x[0, 0] == self.x[-1, 0]:
            for i in range(np.shape(self.mat)[0] - 1):
                for j in range(np.shape(self.mat)[1] - 1):
                    rho[i, j] = self.material_types[self.mat[i, j]].rho
        else:
            for i in range(np.shape(self.mat)[0]):
                for j in range(np.shape(self.mat)[1] - 1):
                    rho[i, j] = self.material_types[self.mat[i, j]].rho
        return rho

    def masses(self):
        dy = self.y[1:] - self.y[:-1]
        volumes = self.Bi[:, :-1] * dy
        rho = self.rho()
        return volumes * rho

    def youngs_mod(self):
        E = np.zeros(np.shape(self.mat))
        for i in range(np.shape(self.mat)[0]):
            for j in range(np.shape(self.mat)[1]):
                E[i, j] = self.material_types[self.mat[i, j]].E
        return E

    def overall_inertia(self):
        dy = self.y[1:] - self.y[:-1]
        masses = self.masses()
        radii = np.sqrt(self.Bi / np.pi)[:, :-1]

        # Apply parallel axis theorem
        x_booms, z_booms = np.split(
            np.reshape(self.section, (np.size(self.y), 2 * np.shape(self.x)[0])), 2, 1
        )
        if np.all(x_booms[:, 0] == x_booms[:, -1]):
            x_booms_nr = x_booms[:, :-1]
            z_booms_nr = z_booms[:, :-1]
        else:
            x_booms_nr = np.copy(x_booms)
            z_booms_nr = np.copy(z_booms)

        x_booms_nr = x_booms_nr.T
        z_booms_nr = z_booms_nr.T

        xcg = np.sum(masses * x_booms_nr[:, :-1], 0) / np.sum(masses, 0)
        zcg = np.sum(masses * z_booms_nr[:, :-1], 0) / np.sum(masses, 0)
        ycg = np.sum(masses * self.y[:-1], 1) / np.sum(masses, 1)

        Ixx = masses * (
            ((self.y[:-1] * np.ones(np.shape(masses))).T - ycg).T ** 2
            + (z_booms_nr[:, :-1] - zcg) ** 2
        )
        Iyy = masses * ((z_booms_nr[:, :-1] - zcg) ** 2 + (x_booms_nr[:, :-1] - xcg) ** 2)
        Izz = masses * (
            ((self.y[:-1] * np.ones(np.shape(masses))).T - ycg).T ** 2
            + (x_booms_nr[:, :-1] - xcg) ** 2
        )

        self.Ix = np.sum(Ixx)
        self.Iy = np.sum(Iyy)
        self.Iz = np.sum(Izz)

        xcg_overall = np.sum(masses * x_booms_nr[:, :-1]) / np.sum(masses)
        zcg_overall = np.sum(masses * z_booms_nr[:, :-1]) / np.sum(masses)
        ycg_overall = np.sum(masses * self.y[:-1]) / np.sum(masses)

        self.xcg = xcg_overall
        self.ycg = ycg_overall
        self.zcg = zcg_overall

        return self.Ix, self.Iy, self.Iz, xcg_overall, ycg_overall, zcg_overall

    def _determine_pitch_and_D(self, n_rows, b):
        # Initial guesses for rivet properties
        D = 0.005
        rivet_mat = materials["Titanium Alloys"]
        n_safety = 3 * 1.5

        # Plate properties
        indx = (np.abs(self.y - b)).argmin()
        width = np.max(self.x[:, indx]) - np.min(self.x[:, indx])
        thickness = np.min(self.t[:, indx])
        plate_mat = self.material_types[self.mat[0, indx]]
        Py = np.abs(self.f_loading[indx][1]) + np.abs(self.m_loading[indx][0] * (np.max(self.section[indx][1]) - self.fix[1][indx]))

        Py = Py / n_rows
        n_rivets = 1
        # Shearing stress - failure of the rivet
        tau = Py / (np.pi * D**2 / 4)
        if tau >= rivet_mat.tau / n_safety:
            D = np.ceil(1000 * np.sqrt(Py / (np.pi * rivet_mat.tau / n_safety / 4))) / 1e3
            if D > 0.02:
                print(
                    f"Rivet diameter is greater than 20 mm - consider increasing plate thickness instead"
                )

        # Bearing stress - failure of the plate
        bearing1 = (Py / n_rivets) / (D * thickness)
        if bearing1 >= plate_mat.compressive / n_safety:
            n_rivets = np.ceil(Py / (plate_mat.compressive / n_safety * D * thickness))

        # Tensions failure - failure of the plate
        w = (width - D * n_rivets) / (n_rivets + 1)
        sigma = Py / (w * thickness)
        while sigma >= plate_mat.tensile / n_safety:
            n_rivets += 1
            w = (width - D * n_rivets) / (n_rivets + 1)
            sigma = Py / (w * thickness)

        return n_rivets, D

    def design_joint(self, b):
        n0, D0 = self._determine_pitch_and_D(1, b)
        n1, D1 = self._determine_pitch_and_D(2, b)

        if n0 * D0 < 2 * n1 * D1:
            print(Fore.BLUE +
                f"The lightest option at {-b} [m] from the root is to have 1 row of {n0} rivets of {D0 * 1e3} mm in diameter"
            )
            print(Fore.WHITE)
            return n0, D0
        else:
            print(Fore.BLUE +
                f"The lightest option at {-b} [m] from the root is to have 2 rows of {n1} rivets of {D1 * 1e3} mm in diameter"
            )
            print(Fore.WHITE)
            return n1, D1

    def _poisson(self):
        v = np.zeros(np.shape((self.Bi[:, :-1] * self.y[1:] - self.y[:-1])))
        if self.x[0, 0] == self.x[-1, 0]:
            for i in range(np.shape(self.mat)[0] - 1):
                for j in range(np.shape(self.mat)[1] - 1):
                    v[i, j] = self.material_types[self.mat[i, j]].poisson
        else:
            for i in range(np.shape(self.mat)[0]):
                for j in range(np.shape(self.mat)[1] - 1):
                    v[i, j] = self.material_types[self.mat[i, j]].poisson
        return v

    def _buckling_stress(self, n_ribs, n_stringers, section):
        E = np.mean(self.youngs_mod(), 0)  # [GPa]
        t = np.mean(self.t, 0)
        v = np.mean(self._poisson(), 0)

        yi = section * int(len(self.y) / (n_ribs + 1))
        yi1 = (section + 1) * int(len(self.y) / (n_ribs + 1))

        b = max(abs(self.x[:, yi])) / (n_stringers + 1)
        Kc = 5.6  # Conservative estimate based on https://core.ac.uk/download/pdf/229099116.pdf

        sig_crit = -Kc * (np.pi**2 / (12 * (1 - np.mean(v[yi:yi1])**2))) * np.mean(E[yi:yi1]) * (np.mean(t[yi:yi1]) / b) ** 2
        return sig_crit


    def buckling(self):
        ## Wingbox Geometry
        y = self.y
        x = self.x

        rho = self.rho()                                        # [kg/m3] The density of the ribs
        A_rib = self.AirfoilArea() * 0.1                        # [m2] Area of the ribs
        t_rib = 0.001                                           # [m] thickness of the ribs
        m_rib = np.hstack((rho[0], rho[0, 0])) * A_rib * t_rib  # [kg] mass of a rib
        A_str = 0.04 * 0.001                                    # [m2] Area of the stringers

        sigma = self.sigma
        iteration_mass = []
        mass = 0
        ItInfo = []
        if np.any(sigma < 0):
            for i in tqdm(range(int(len(y) / 2))):  # Guess number of ribs
                N_ribs = i
                strdis = []
                sig = []
                dy = int(len(y) / (N_ribs + 1)) - 1
                for j in range(i + 1):  # Iterate over sections
                    sig_crit = 1e10  # Initial guess - high so that it always goes into the while loop
                    k = -1
                    yi = j * int(len(y) / (N_ribs + 1))
                    yi1 = (j + 1) * int(len(y) / (N_ribs + 1))
                    while np.any(sigma[:, yi:yi1] < sig_crit):
                        k += 1
                        sig_crit = self._buckling_stress(n_ribs=N_ribs, n_stringers=k, section=j)
                    if k == -1:
                        k = 0
                    mass += (
                            k * A_str * np.min(rho) * abs(y[int(j * dy)] - y[int(((j + 1) * dy))])
                            + m_rib[int(((j + 1) * dy))]
                    )
                    strdis.append(k)
                    sig.append(sig_crit)
                iteration_mass.append(mass)
                ItInfo.append([mass, strdis, sig])
            np.array(iteration_mass)
            I = np.where(iteration_mass[1:] == np.min(iteration_mass[1:]))[0][0]
            print(Fore.BLUE +
                  f"The best stringer distribution is {ItInfo[I + 1][1]}, having a mass of {ItInfo[I + 1][0]} [kg]"
                  )
            return ItInfo[I + 1]
        else:
            print(Fore.BLUE + 'All members are in tension. No buckling will occur' + Fore.WHITE)
            return ItInfo

    def plot_internal_loading(self):
        fig, (axs1, axs2) = plt.subplots(2, 1)
        f_loading = np.reshape(self.f_loading, (len(self.y), 3))
        m_loading = np.reshape(self.m_loading, (len(self.y), 3))

        axs1.plot(self.y, f_loading[:, 0] / 1e3, label=r"$V_x$")
        axs1.plot(self.y, f_loading[:, 1] / 1e3, label=r"$V_y$")
        axs1.plot(self.y, f_loading[:, 2] / 1e3, label=r"$V_z$")
        axs1.set_xlabel("Span [m]")
        axs1.set_ylabel("Force [kN]")
        axs1.legend()

        axs2.plot(self.y, m_loading[:, 0] / 1e3, label=r"$M_x$")
        axs2.plot(self.y, m_loading[:, 1] / 1e3, label=r"$M_y$")
        axs2.plot(self.y, m_loading[:, 2] / 1e3, label=r"$M_z$")
        axs2.set_xlabel("Span [m]")
        axs2.set_ylabel("Moment [kNm]")
        axs2.legend()

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    ...
