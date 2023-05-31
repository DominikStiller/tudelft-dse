from dse.detailed.Structures.material_properties import materials
from scipy.interpolate import InterpolatedUnivariateSpline
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
    indx = Data.index(["Main Wing"])
    Data = Data[indx + 1 :]

    xflr_data = pd.DataFrame(Data[1:], columns=Data[0], dtype=float)

    # Cut half wing
    cutpoint = int(round(len(xflr_data["Cl"]) / 2))
    xflr_data = xflr_data[:cutpoint]

    ## Application of the forces
    # Assume drag acts through the neutral axis
    application = np.zeros((3, len(xflr_data['  y-span'])))   # [m] (3, n) where n is the length of y-span
    application[0] = np.array(xflr_data['XCP'])             # [m] Point of application in the x direction
    application[1] = np.array(xflr_data['  y-span'])          # [m] Point of application in the y direction

    ## Magnitude of the forces
    # Overwrite the forces in y direction if there are any present.
    # [m2] Array of all the surface areas
    S_array = xflr_data["Chord"] * np.hstack(
        (
            np.array([b - abs(xflr_data["  y-span"][0])]),
            abs(application[1][1:] - application[1][:-1]),
        )
    )
    mag = np.zeros((3, len(xflr_data["  y-span"])))  # [N] (3 x n) where n is the length of y_lst
    mag[0] = (
        -(xflr_data["ICd"] + xflr_data["PCd"]) * q * S_array
    )  # [N] Forces in x-direction due to drag
    mag[2] = xflr_data["Cl"] * q * S_array  # [N] Forces in z-direction due to lift

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
    def __init__(self, width, length, height, cross_section, material):
        # Beam dimensions
        if type(length) == int or type(length) == float:
            self.y = np.linspace(-length, 0, 100)
        elif type(length) == np.ndarray:
            self.y = length
        else:
            raise TypeError("Length needs to be either a constant float/int or a 1D array")

        if type(width) == int or type(width) == float:
            self.x = np.reshape(np.linspace(-width/2, width, 100), (100, 1)) * np.size(self.y)
        elif type(width) == np.ndarray:
            self.x = width
        else:
            raise TypeError("Width needs to be either a constant float/int or a 100 x n array")


        if type(height) == int or type(height) == float:
            self.z = np.reshape(np.linspace(-height/2, height/2, 100), (100, 1)) * np.ones(np.size(self.y))
        elif type(height) == np.ndarray:
            self.z = height
        else:
            raise TypeError("Height needs to be either a constant float/int or a 100 x n array")

        if np.shape(self.x) != np.shape(self.z):
            raise ValueError('Width and height must have matching shapes')

        if type(cross_section) == str:
            if cross_section == "full":
                self.section = np.ones((len(self.y), len(self.z), len(self.x)))
            elif cross_section == "constant":
                self.section = np.ones((len(self.y), 2, len(self.x))) * np.vstack((self.x, self.z))
        else:
            if type(cross_section) != np.ndarray:
                raise TypeError('The cross section needs to be either "full" or a 3D array')
            else:
                if np.shape(cross_section) != (len(self.y), 2, len(self.x)):
                    raise ValueError('Cross section needs to have a size consistent with the '
                                     'amount of points indicated for length, height and width')
                if len(np.shape(self.x)) != 2:
                    raise ValueError('If cross_section != "full" or "constant", then width and height need to be 2D '
                                     'arrays')
                self.section = cross_section

        # Internal forces
        self.f_loading = np.zeros((len(self.y), 3, 1))

        # Internal moments
        self.m_loading = np.zeros((len(self.y), 3, 1))

        self.mat = material

        # Parameters to be calulated later
        self.t = None
        self.sigma = None

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
                    [[np.mean(self.x), self.y[y_point + j], np.mean(self.z)]]
                ).T - np.reshape(force.application[:, i], (3, 1))
                self.m_loading[y_point + j] += np.reshape(
                    np.cross(force.F[:, i], distance.T), (3, 1)
                )

    def AirfoilBoom(self):
        ## The initial guess of the skin thickness follows from the inful percentage of the wingbox, and  circumference
        Infill = 0.10  # The percentage infull of the wingbox area
        # Locating the leading edge.
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

        ## Reshapes the x and z coordinates to the correct shapes, one with repeating first node, and one without (_nr)
        x_booms = np.reshape(self.x, (len(self.x), 1))
        z_booms = np.reshape(self.z, (len(self.z), 1))
        x_booms_nr = x_booms[:-1]
        z_booms_nr = z_booms[:-1]

        Bi_initial = A_airfoil / len(x_booms_nr) * Infill

        return Bi_initial

    def unitSquareBoom(self):
        ### Calculates the initial guess for the boom areas for a square wingbox of unit length
        ## The initial guess of the skin thickness follows from the inful percentage of the wingbox, and  circumference

        Infill = 0.10  # The percentage infull of the wingbox area
        A_airfoil = 1

        ## Reshapes the x and z coordinates to the correct shapes, one with repeating first node, and one without (_nr)
        x_booms = np.reshape(self.x, (len(self.x), 1))
        z_booms = np.reshape(self.z, (len(self.z), 1))
        x_booms_nr = x_booms[:-1]
        z_booms_nr = z_booms[:-1]

        Bi_initial = A_airfoil / len(x_booms_nr) * Infill

        return Bi_initial

    def DesignConstraints(self):
        n = 1.5  # [-] Safety factor
        sigma_yield = self.mat.sigmay  # [MPa] The yield strength
        tSkin_min = 0.001  # [m] The minimal allowable thickness
        return n, sigma_yield, tSkin_min

    def NeutralAxis(self, boomArea_nr, x_booms_nr, z_booms_nr):
        # Neutral Axis Calculations
        NAx = (np.sum(boomArea_nr * x_booms_nr, 0)) / (np.sum(boomArea_nr, 0))  # [m]
        NAz = (np.sum(boomArea_nr * z_booms_nr, 0)) / (np.sum(boomArea_nr, 0))  # [m]
        return NAx, NAz

    def MoI(self, boomArea_nr, x_booms_nr, z_booms_nr):
        # Moment of Inertia calculations
        NAx, NAz = self.NeutralAxis(boomArea_nr, x_booms_nr, z_booms_nr)
        Ixx = np.sum(boomArea_nr * (z_booms_nr - NAz) ** 2, 0)
        Izz = np.sum(boomArea_nr * (x_booms_nr - NAx) ** 2, 0)
        Ixz = np.sum(boomArea_nr * (x_booms_nr - NAx) * (z_booms_nr - NAz), 0)
        return Ixx, Izz, Ixz

    def StressCalculations(self, boomArea_nr):
        Mx = np.reshape(self.m_loading[:, 0], np.size(self.y))
        Mz = np.reshape(self.m_loading[:, 2], np.size(self.y))
        Fy = np.reshape(self.f_loading[:, 1], np.size(self.y))

        x_booms, z_booms = np.split(np.reshape(self.section, (np.size(self.y), 2 * np.shape(self.x)[0])), 2, 1)
        if np.all(x_booms[:, 0] == x_booms[:, -1]):
            x_booms = x_booms[:, :-1]
            z_booms = z_booms[:, :-1]

        x_booms = x_booms.T
        z_booms = z_booms.T

        # Call the neutral & moments of inertia
        NAx, NAz = self.NeutralAxis(boomArea_nr, x_booms, z_booms)
        Ixx, Izz, Ixz = self.MoI(boomArea_nr, x_booms, z_booms)
        # Stress calculations
        sigma_nr = (
                           (Mx * Izz - Mz * Ixz) * (z_booms - NAz) + (Mz * Ixx - Mx * Ixz) * (x_booms - NAx)
                   ) / (Ixx * Izz + Ixz**2) + (Fy / np.shape(boomArea_nr)[0]) / (boomArea_nr)
        return sigma_nr

    def BoomArea(self, boomAreaCopy, tSkin, boomDistance, sigma):
        # Update the boom areas for the given stress.
        for i in range(np.shape(boomAreaCopy)[0]):
            boomAreaCopy[i] = tSkin[i - 1] * boomDistance[i - 1] / 6 * (
                    2 + sigma[i - 1] / sigma[i]
            ) + tSkin[i] * boomDistance[i] / 6 * (2 + sigma[i + 1] / sigma[i])
        return boomAreaCopy

    def InternalStress(self, boom_coordinates, interconnection, i):
        if interconnection != 0:  # define the interconnection of all of the boom areas
            print("Interconnection between stringers still need to be implemented")

        if boom_coordinates != 0:  # Defines locations of the booms
            print("The code now only runs for all coordinates of the stringers being the booms")

        x_booms, z_booms = np.split(np.reshape(self.section, (np.size(self.y), 2 * np.shape(self.x)[0])), 2, 1)
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
        if i == 0:
            Bi_initial = self.AirfoilBoom()
        elif i == 1:
            Bi_initial = self.unitSquareBoom()
        else:
            Bi_initial = 0.001 * np.ones(np.shape(x_booms_nr))
        boomArea_nr = np.ones((np.shape(x_booms_nr)[0], np.size(self.y))) * Bi_initial
        tSkin = np.ones(np.shape(boomDistance)) * tSkin_min
        sigma = np.ones(np.shape(self.x)) * n * sigma_ult

        while np.any(np.abs(sigma - sigma_ult / n) / (sigma_ult / n) > 0.1):
            diff = np.ones(np.shape(boomArea_nr)[1])
            while np.any(np.abs(diff) > 0.0001):
                # Stress Calculations
                sigma_nr = self.StressCalculations(boomArea_nr)

                # Stress with repeating column for the first node for easy calculations for the boom areas
                sigma = np.vstack((sigma_nr, sigma_nr[0]))
                sigma[np.where(np.abs(sigma) <= 0.01)] = 0.1

                # Boom area calculations
                boomAreaCopy = np.copy(boomArea_nr)
                boomAreaCopy = self.BoomArea(boomAreaCopy, tSkin, boomDistance, sigma)

                # Calculating the difference between the new and old boom areas.
                diff = np.mean((boomAreaCopy - boomArea_nr) / boomArea_nr, 0)
                boomArea_nr = np.copy(boomAreaCopy)

            indx = 0
            for k in range(np.shape(sigma)[1]):
                if np.max(abs(sigma[:, k])) >= sigma_ult:
                    indx += 1
            if indx != 0:
                tSkin[:, -indx:] += 0.001
            else:
                break
            if np.any(np.isnan(sigma)):
                raise ValueError('There was a NaN in the stress array')

        self.t = tSkin
        self.sigma = sigma_nr
        # plt.axhline(sigma_ult)
        # plt.axhline(-sigma_ult)
        plt.plot(self.y, sigma.transpose())
        plt.show()

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
    g = 3.71
    q = 0.5 * 0.01 * 112 ** 2
    aerodynamic_forces = xflr_forces('Test_xflr5_file.csv', q, 16.8)

    chord = 3.35
    Airfoil = pd.read_csv("S1223.dat", delimiter="\s+", dtype=float, skiprows=1, names=["x", "z"])
    l = np.linspace(-16.8, 0, 100)
    wing = Beam(
        width=Airfoil["x"].to_numpy() * chord,
        height=Airfoil["z"].to_numpy() * chord,
        length=l,
        cross_section="constant",
        material=materials['Al/Si']
    )
    # wing.add_loading(aerodynamic_forces)
    # thrust = Force(
    #     magnitude=np.array(
    #         [[739.5/2],
    #          [0],
    #          [0]]
    #     ),
    #     point_of_application=np.array(
    #         [[0],
    #          [-16.8],
    #          [0]]
    #     )
    # )
    fuselage_height = 2
    theta = np.arctan(fuselage_height/16.8)

    MTOM = 3000
    m_r = 910 / 2
    m_e = 100 / 2

    bracing_TO = Force(
        magnitude=np.array(
            [
                [0],
                [g / np.tan(theta) * (1.1 * MTOM / 2 - (m_r + m_e))],
                [-g * (1.1 * MTOM / 2 - (m_r + m_e))]
            ]
        ),
        point_of_application=np.array(
            [
                [1.603],
                [-16.8],
                [0.1742]
            ]
        )
    )
    engine_and_rotor_weight = Force(
        magnitude=np.array(
            [
                [0],
                [0],
                [-(m_r + m_e) * g]
            ]
        ),
        point_of_application=np.array(
            [
                [1.603],
                [-16.8],
                [0.1742]
            ]
        )
    )
    liftOffLoad = Force(
        magnitude=np.array(
            [
                [0],
                [0],
                [1.1 * MTOM * g / 2]
            ]
        ),
        point_of_application=np.array(
            [
                [1.603],
                [-16.8],
                [0.1742]
            ]
        )
    )
    wing.add_loading(liftOffLoad)
    wing.plot_internal_loading()

    wing.add_loading(engine_and_rotor_weight)
    wing.plot_internal_loading()

    wing.add_loading(bracing_TO)
    wing.plot_internal_loading()

    wing.InternalStress(0, 0, 0)

    plt.plot(wing.y, np.max(wing.t, 0) * 1000, label='Max thickness')
    plt.plot(wing.y, np.min(wing.t, 0) * 1000, label='Min thickness')
    plt.xlabel('Span [m]')
    plt.ylabel('Thickness [mm]')
    plt.legend()
    plt.show()
