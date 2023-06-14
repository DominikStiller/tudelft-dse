from numpy.testing import assert_allclose
from dse.detailed.Structures.StructureClasses import *
from dse.detailed.Structures.material_properties import materials
from unittest import TestCase


class TestBeam(TestCase):
    # Functions that need testing:
    # TorsionStress
    # BoomArea
    # Buckling
    def test_add_loading_point_load(self):
        x = 0.1
        y = 2
        z = 0.05
        wing = Beam(x, y, z, "square", material='Al/Si',
                    fixing_points=np.array([[x/2], [z/2]]) * np.ones(100))

        F = 1
        loc = -1.5

        point_load = Force(
            magnitude=np.array([[0, 0, F]]).T, point_of_application=np.array([[x/2, loc, z/2]]).T
        )

        wing.add_loading(point_load)

        y_index = (np.abs(wing.y - loc)).argmin()
        assert_allclose(((wing.m_loading[y_index] - wing.m_loading[-1]) / -loc)[0], F, rtol=0.01)

    def test_add_loading_distributed_load(self):
        x = 0.1
        y = 2
        z = 0.05
        wing = Beam(x, y, z, "square", material='Al/Si',
                    fixing_points=np.array([[x/2], [z/2]]) * np.ones(100))

        F = np.array([[0, 0, 0], [0, 0, 0], [1, 1, 1]])

        app = np.array([[x/2, x/2, x/2], [-0.1, -1, -1.5], [z/2, z/2, z/2]])

        force = Force(magnitude=F, point_of_application=app)
        wing.add_loading(force)

        y_index_1 = (np.abs(wing.y - app[:, 0][1])).argmin()
        y_index_2 = (np.abs(wing.y - app[:, 1][1])).argmin()
        y_index_3 = (np.abs(wing.y - app[:, 2][1])).argmin()

        # Check moments
        assert_allclose(wing.m_loading[-1][0], -2.6, rtol=0.05, err_msg='Moment at the tip is wrong')
        assert_allclose(wing.m_loading[y_index_1][0], -2.3, rtol=0.05, err_msg='Moment at segment 1 is wrong')
        assert_allclose(wing.m_loading[y_index_2][0], -0.5, rtol=0.05, err_msg='Moment at segment 2 is wrong')
        assert_allclose(wing.m_loading[y_index_3-1][0], 0, rtol=0.05, err_msg='Moment at the third point is wrong')

        # Check forces
        self.assertTrue(np.all(wing.f_loading[y_index_1:] == np.reshape(np.sum(F, 1), (3, 1))))
        self.assertTrue(np.all(wing.f_loading[y_index_2:y_index_1] == np.reshape(np.sum(F[:, 1:], 1), (3, 1))))
        self.assertTrue(np.all(wing.f_loading[y_index_3:y_index_2] == np.reshape(F[:, -1], (3, 1))))

    def test_InternalLoads_Rhombus(self):
        ## Hand Calculations
        NAx_hand = 0.  # Neutral axis as obtained by hand calculations
        NAz_hand = 0.  # Neutral axis as obtained by hand calculations
        Ixx_hand = 0.0375  # [m4]
        Izz_hand = 0.0375  # [m4]
        Ixz_hand = 0.000   # [m4]
        sigma_hand = np.array(
            [
                [0.266666666666, 1.4666666666, 2.666666666, 1.2, -0.266666666666, -1.4666666666, -2.6666666, -1.2]
            ]
        ).T * 1e6
        b = 10
        x = np.array([[1, 0.5, 0, -0.5, -1, -0.5, 0, 0.5]])
        z = np.array([[0, 0.5, 1, 0.5, 0, -0.5, -1, -0.5]])
        l = np.linspace(-b, 0, 100)
        rombus = Beam(
            width=x.T * np.ones(100),
            height=z.T * np.ones(100),
            length=l,
            cross_section=np.vstack((x, z)) * np.ones((np.size(l), 2, 1)),
            material='Al/Si',
            fixing_points=np.array([[0], [0]]) * np.ones(100)
        )
        Applied_Load = Force(
            magnitude=np.array(
                [
                    [1000],
                    [0],
                    [-10000]
                ]
            ),
            point_of_application=np.array(
                [
                    [0],
                    [-b],
                    [0]
                ]
            )
        )
        rombus.add_loading(Applied_Load)
        # rombus.plot_internal_loading()

        x_booms_nr, z_booms_nr = np.split(np.reshape(rombus.section, (np.size(rombus.y), 2 * np.shape(rombus.x)[0])), 2,
                                          1)
        if np.all(x_booms_nr[:, 0] == x_booms_nr[:, -1]):
            x_booms_nr = x_booms_nr[:, :-1]
            z_booms_nr = x_booms_nr[:, :-1]

        x_booms_nr = x_booms_nr.T
        z_booms_nr = z_booms_nr.T

        NAx, NAz = rombus.NeutralAxis(boomArea_nr=np.ones((8, np.size(l)))*0.0125, x_booms_nr=x_booms_nr, z_booms_nr=z_booms_nr)
        Ixx, Izz, Ixz = rombus.MoI(boomArea_nr=np.ones((8, np.size(l)))*0.0125, x_booms_nr=x_booms_nr, z_booms_nr=z_booms_nr)
        sigma_nr = rombus.StressCalculations(boomArea_nr=np.ones((8, np.size(l)))*0.0125)

        assert_allclose(NAx[-1], NAx_hand, atol=1e-5, err_msg='The NAx is off')
        assert_allclose(NAz[-1], NAz_hand, atol=1e-5, err_msg='The NAz is off')
        assert_allclose(Ixx[-1], Ixx_hand, rtol=1e-4, err_msg='The Ixx is off')
        assert_allclose(Izz[-1], Izz_hand, rtol=1e-4, err_msg='The Izz is off')
        assert_allclose(Ixz[-1], Ixz_hand, atol=1e-5, err_msg='The Ixz is off')
        assert_allclose(sigma_nr[:, -1], sigma_hand.T[0], rtol=1e-3, err_msg='The stresses are wrong')

    def test_InternalLoads_Parallelogram(self):
        ## Hand Calculations
        NAx_hand = 0.  # Neutral axis as obtained by hand calculations
        NAz_hand = 0.  # Neutral axis as obtained by hand calculations
        Ixx_hand = 0.0375  # [m4]
        Izz_hand = 0.075  # [m4]
        Ixz_hand = 0.0375  # [m4]
        sigma_hand = np.array(
            [
                [-2.4, 0.133333, 2.6666, 2.5333, 2.4, -0.13333, -2.666, -2.53333]
            ]
        ).T * 1e6
        b = 10
        x = np.array([[1, 1, 1, 0, -1, -1, -1, 0]])
        z = np.array([[0, 0.5, 1, 0.5, 0, -0.5, -1, -0.5]])
        l = np.linspace(-b, 0, 100)
        parallelogram = Beam(
            width=x.T * np.ones(np.size(l)),
            height=z.T * np.ones(np.size(l)),
            length=l,
            cross_section=np.vstack((x, z)) * np.ones((np.size(l), 2, 1)),
            material='Al/Si',
            fixing_points=np.array([[0], [0]]) * np.ones(np.size(l))
        )
        Applied_Load = Force(
            magnitude=np.array(
                [
                    [1000],
                    [0],
                    [-10000]
                ]
            ),
            point_of_application=np.array(
                [
                    [0],
                    [-b],
                    [0]
                ]
            )
        )
        parallelogram.add_loading(Applied_Load)

        x_booms_nr, z_booms_nr = np.split(np.reshape(parallelogram.section, (np.size(parallelogram.y), 2 * np.shape(parallelogram.x)[0])), 2,
                                          1)
        if np.all(x_booms_nr[:, 0] == x_booms_nr[:, -1]):
            x_booms_nr = x_booms_nr[:, :-1]
            z_booms_nr = x_booms_nr[:, :-1]

        x_booms_nr = x_booms_nr.T
        z_booms_nr = z_booms_nr.T

        NAx, NAz = parallelogram.NeutralAxis(boomArea_nr=np.ones((8, np.size(l)))*0.0125, x_booms_nr=x_booms_nr, z_booms_nr=z_booms_nr)
        Ixx, Izz, Ixz = parallelogram.MoI(boomArea_nr=np.ones((8, np.size(l)))*0.0125, x_booms_nr=x_booms_nr, z_booms_nr=z_booms_nr)
        sigma_nr = parallelogram.StressCalculations(boomArea_nr=np.ones((8, np.size(l)))*0.0125)

        assert_allclose(NAx[-1], NAx_hand, atol=1e-5, err_msg='The NAx is off')
        assert_allclose(NAz[-1], NAz_hand, atol=1e-5, err_msg='The NAz is off')
        assert_allclose(Ixx[-1], Ixx_hand, rtol=1e-4, err_msg='The Ixx is off')
        assert_allclose(Izz[-1], Izz_hand, rtol=1e-4, err_msg='The Izz is off')
        assert_allclose(Ixz[-1], Ixz_hand, atol=1e-5, err_msg='The Ixz is off')
        assert_allclose(sigma_nr[:, -1], sigma_hand.T[0], rtol=1e-3, err_msg='The stresses are wrong')

    def test_ShearLoads_Wingbox(self):
        ## Hand Calculations
        q_total_hand = np.array(
            [
                [12.72, -5.32, -34.18, -37.79, -34.18, -5.32, 12.72, 17.05]
            ]
        ).T * 1e3 * np.ones(100)
        b = 10
        x = np.array([[0.6, 0.36, 0.12, 0, 0, .12, .36, .6, .6]])
        z = np.array([[0.03, 0.1, .1, 0.05, -0.05, -0.1, -0.1, -0.03, 0.03]])
        l = np.linspace(-b, 0, 100)
        boomArea_nr = np.array([[0.0002, 0.00025, 0.0004, 0.0001, 0.0001, 0.0004, 0.00025, 0.0002,]]).T * np.ones(len(l))
        wingbox = Beam(
            width=x.T * np.ones(np.size(l)),
            height=z.T * np.ones(np.size(l)),
            length=l,
            cross_section=np.vstack((x, z)) * np.ones((np.size(l), 2, 1)),
            material='Al/Si',
            fixing_points=np.array([[0], [0]]) * np.ones(np.size(l))
        )
        Applied_Load = Force(
            magnitude=np.array(
                [
                    [0],
                    [0],
                    [10000]
                ]
            ),
            point_of_application=np.array(
                [
                    [0.12],
                    [-b],
                    [0]
                ]
            )
        )
        wingbox.add_loading(Applied_Load)

        x_booms_nr, z_booms_nr = np.split(
            np.reshape(wingbox.section, (np.size(wingbox.y), 2 * np.shape(wingbox.x)[0])), 2,
            1)
        if np.all(x_booms_nr[:, 0] == x_booms_nr[:, -1]):
            x_booms_nr = x_booms_nr[:, :-1]
            z_booms_nr = x_booms_nr[:, :-1]


        q_total = wingbox.TorsionStress(boomArea_nr=boomArea_nr, A_Airfoil=97200 * 1e-6)
        assert_allclose(q_total, q_total_hand, rtol=1e-3, err_msg='The stresses are wrong')

    def test_stress_calculations(self):
        x = np.atleast_2d(np.hstack(
            (np.linspace(-1, 0, 25)[:-1],
             np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25)[:-1])
        ))
        z = np.atleast_2d(np.hstack(
            (np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25)[:-1],
             np.linspace(-1, 0, 25)[:-1]
             )
        ))
        rombus = Beam(
            width=x.T * np.ones(100),
            height=z.T * np.ones(100),
            length=5,
            cross_section=np.vstack((x, z)) * np.ones((100, 1, 1)),
            material='Al/Si',
            fixing_points=np.array([[0], [0]]) * np.ones(100)
        )
        A = np.ones((np.size(x), 100))

        dummy_force_1 = Force(
            magnitude=np.array(
                [
                    [0],
                    [1e3],
                    [0]
                ]
            ),
            point_of_application=np.array(
                [
                    [0],
                    [-5],
                    [0]
                ]
            )
        )
        rombus.add_loading(dummy_force_1)

        stress = rombus.StressCalculations(
            boomArea_nr=A
        )

        self.assertTrue(np.all(stress == stress[0]), msg='The stress on the booms due to a point load along the length '
                                                         'does not result in a constant stress')
        rombus.unload()

        dummy_force_2 = Force(
            magnitude=np.array(
                [
                    [0],
                    [0],
                    [1e3]
                ]
            ),
            point_of_application=np.array(
                [
                    [0],
                    [-5],
                    [0]
                ]
            )
        )

        rombus.add_loading(dummy_force_2)

        stress = rombus.StressCalculations(
            boomArea_nr=A
        )

        stress_split = np.split(stress, 4)

        assert_allclose(stress_split[0][1:], np.flip(stress_split[1], 0)[:-1], rtol=1e-4,
                        err_msg='The top booms have non-symmetrical stresses')
        assert_allclose(stress_split[2][1:], np.flip(stress_split[3], 0)[:-1], rtol=1e-4,
                        err_msg='The bottom booms have non-symmetrical stresses')
        self.assertTrue(np.where(stress == np.min(stress))[0] == 24,
                        msg='The maximum compressive stress is not at the top')
        self.assertTrue(np.where(stress == np.max(stress))[0] == 72,
                        msg='The maximum tensile stress is not at the bottom')
        assert_allclose(np.max(stress), -np.min(stress), rtol=1e-4,
                        err_msg='The maximum tensile and compressive stresses are not the same')

    def test_internal_stress(self):
        x = np.atleast_2d(np.hstack(
            (np.linspace(-1, 0, 25)[:-1],
             np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25))
        ))
        z = np.atleast_2d(np.hstack(
            (np.linspace(0, 1, 25)[:-1],
             np.linspace(1, 0, 25)[:-1],
             np.linspace(0, -1, 25)[:-1],
             np.linspace(-1, 0, 25)
             )
        ))
        rombus = Beam(
            width=x.T * np.ones(100),
            height=z.T * np.ones(100),
            length=5,
            cross_section=np.vstack((x, z)) * np.ones((100, 1, 1)),
            material='Al/Si',
            fixing_points=np.array([[0], [0]]) * np.ones(100)
        )

        dummy_force_2 = Force(
            magnitude=np.array(
                [
                    [0],
                    [0],
                    [1e3]
                ]
            ),
            point_of_application=np.array(
                [
                    [0],
                    [-5],
                    [0]
                ]
            )
        )

        A = (np.max(x) - np.min(x)) * (np.max(z) - np.min(z)) / 2

        rombus.add_loading(dummy_force_2)
        rombus.InternalStress(0, 0, A)

        self.assertTupleEqual(np.shape(rombus.t), (np.size(x)-1, 100),
                              msg='rombus.t is not an array of the correct shape')
        self.assertTrue(np.all(rombus.t >= 0.001), msg='The minimum thickness is smaller than 1mm')
        self.assertTrue(np.all(rombus.sigma <= 250e6/1.5), msg='Maximum stress is higher than allowed')

    def test_masses_constant_section(self):
        x = 1
        z = 1
        y = 2

        squareBeam = Beam(
            width=x,
            height=z,
            length=y,
            cross_section='square',
            material='Al/Si',
            fixing_points=np.array([[x/2], [z/2]]) * np.ones(100)
        )

        infill = 0.1
        A = x * z * infill
        squareBeam.Bi = A / np.shape(squareBeam.x)[0] * np.ones(np.shape(squareBeam.x))
        masses = squareBeam.masses()

        self.assertTrue(np.all(masses == masses[0]), msg='Discrete masses are not constant')
        assert_allclose(np.sum(masses), A * y * materials['Al/Si'].rho, rtol=1e-4, err_msg='Total mass is wrong')

    def test_masses_linear_section(self):
        x = np.atleast_2d(np.hstack((
                    np.zeros(20),
                    np.linspace(0, 1, 21)[:-1],
                    np.ones(20)*1,
                    np.linspace(1, 0, 21)[:-1]
                )))
        z = np.atleast_2d(np.hstack((
                    np.linspace(0, 1, 21)[:-1],
                    np.ones(20) * 1,
                    np.linspace(1, 0, 21)[:-1],
                    np.zeros(20)
                )))
        y = np.linspace(-2, 0, 100)
        piramid = Beam(
            width=x.T * np.ones(100),
            height=z.T * np.linspace(1, 2, 100),
            length=y,
            cross_section=np.vstack((x, z)) * np.ones((100, 1, 1)),
            material='Al/Si',
            fixing_points=np.array([[1], [1]]) * np.ones(100)
        )
        infill = 0.1
        A = np.max(piramid.x, 0) * np.max(piramid.z, 0) * infill
        piramid.Bi = A / np.shape(piramid.x)[0] * np.ones(np.shape(piramid.x))
        masses = piramid.masses()

        masses_y = np.sum(masses, 0)
        summed_masses = np.zeros(np.shape(masses_y))
        for i in range(len(masses_y)-1):
            summed_masses[i] = np.sum(masses_y[:i+1])
        hand_masses_y = materials['Al/Si'].rho * np.flip(np.abs(y) + 0.25 * y**2) * infill

        self.assertTrue(np.all(masses[:, 1:]/masses[:, :-1] > 1), msg='Mass is not increasing')
        assert_allclose(summed_masses[:-1], hand_masses_y[1:-1], rtol=0.01,
                        err_msg="Mass array doesn't coincide with the analytical calculation")

    def test_overall_inertia(self):
        x0 = 1
        z0 = 1
        y = 5

        squareBeam = Beam(
            width=x0,
            height=z0,
            length=y,
            cross_section='square',
            material='Al/Si',
            fixing_points=np.array([[x0/2], [z0/2]]) * np.ones(100)
        )
        t = 0.005
        x1 = x0 - 2*t
        z1 = z0 - 2*t

        A = x0*z0 - x1*z1
        squareBeam.Bi = A / np.shape(squareBeam.x)[0] * np.ones(np.shape(squareBeam.x))
        squareBeam.overall_inertia()

        # Hand calculations
        M0 = x0*z0*y*materials['Al/Si'].rho
        M1 = x1*z1*y*materials['Al/Si'].rho
        Ixx = M0*(y**2 + z0**2)/12 - M1*(y**2 + z1**2)/12
        Izz = M0*(y**2 + x0**2)/12 - M1*(y**2 + x1**2)/12
        Iyy = M0*(x0**2 + z0**2)/12 - M1*(x1**2 + z1**2)/12

        assert_allclose(squareBeam.Ix, Ixx, rtol=1e-3, err_msg='Ixx is wrong')
        assert_allclose(squareBeam.Iy, Iyy, rtol=0.015, err_msg='Iyy is wrong')
        assert_allclose(squareBeam.Iz, Izz, rtol=1e-3, err_msg='Izz is wrong')

    def test_design_joint(self):
        x0 = 1
        z0 = 1
        y = 5

        squareBeam = Beam(
            width=x0,
            height=z0,
            length=y,
            cross_section='square',
            material='Al/Si',
            fixing_points=np.array([[x0 / 2], [z0 / 2]]) * np.ones(100)
        )

        t_plate = 0.001
        squareBeam.t = t_plate * np.ones(np.shape(squareBeam.x))
        P_mag = 1e5
        P = Force(magnitude=np.array([[0, 0, P_mag]]).T,
                  point_of_application=np.array([[x0 / 2], [-y], [z0 / 2]]))

        squareBeam.add_loading(P)
        indx = (np.abs(squareBeam.y + y/2)).argmin()
        F = np.abs(squareBeam.f_loading[indx][1]) + np.abs(squareBeam.m_loading[indx][0]) * np.max(np.abs(squareBeam.z[:, indx] - squareBeam.fix[1, indx]))


        tau_max = materials['Titanium Alloys'].tau / 4.5
        D = np.ceil(1e3*np.sqrt(4 * F / (np.pi * tau_max))) / 1e3

        sigma_max = materials['Al/Si'].compressive / 4.5
        n0 = np.ceil(F / (sigma_max * t_plate * D))

        n1 = n0
        width = np.max(squareBeam.x[:, indx]) - np.min(squareBeam.x[:, indx])
        w = (width - D * n1) / (n1 + 1)
        while F / (t_plate * w) > sigma_max:
            n1 += 1
            w = (width - D * n1) / (n1 + 1)

        n_out, D_out = squareBeam.design_joint(-y/2)

        assert n_out == n1, 'The number of rivets does not coincide with the analytical solution'
        assert D_out == D, 'Rivet diameter does not coincide with the analytical solution'


class ValidateBeam(TestCase):
    def test_internal_load(self):
        # From https://www.scribd.com/document/263390031/Experiment-Shearing-Forces#
        loads0 = np.array([
            [-10],
            [-10],
            [-10]
        ])

        loads1 = np.array([
            [-5, -5, -5],
            [-10, -10, -10],
            [-15, -15, -15],
        ])

        positions0 = np.array([
            [-0.6],
            [-0.3],
            [-0.15]
        ])

        positions1 = np.array([
            [-0.15, -0.3, -0.6],
            [-0.15, -0.3, -0.6],
            [-0.15, -0.3, -0.6]
        ])

        Ra0 = np.array([
            [3.33],
            [6.67],
            [8.33]
        ])

        Ra1 = np.array([
            [9.17],
            [18.33],
            [27.5]
        ])

        Rb0 = np.array([[-6.67, -3.33, -1.67]]).T
        Rb1 = np.array([[-5.83, -11.67, -17.5]]).T

        beam = Beam(
            width=0.05,
            height=0.05,
            length=0.9,
            cross_section='square',
            material='Al/Si',
            fixing_points=np.array([[0.025], [0.025]])*np.ones(100)
        )

        for i in range(3):
            beam.unload()
            loading0 = Force(
                magnitude=np.vstack((np.zeros((2, 1)), loads0[i][0])),
                point_of_application=np.array([
                    [0.025],
                    [positions0[i][0]],
                    [0.025]
                ])
            )
            reaction = Force(
                magnitude=np.array([[0, 0, Ra0[i][0]]]).T,
                point_of_application=np.array([[0.025, -0.8999, 0.025]]).T
            )
            beam.add_loading(loading0)
            beam.add_loading(reaction)
            assert_allclose(beam.f_loading[-1][-1], Rb0[i], rtol=0.01, err_msg='Reaction force does not match')

            beam.unload()
            loading1 = Force(
                magnitude=np.vstack((np.zeros((2, 3)), loads1[i])),
                point_of_application=np.vstack((
                    0.025*np.ones(3),
                    positions1[i],
                    0.025*np.ones(3)
                ))
            )
            reaction = Force(
                magnitude=np.array([[0, 0, Ra1[i][0]]]).T,
                point_of_application=np.array([[0.025, -0.8999, 0.025]]).T
            )
            beam.add_loading(loading1)
            beam.add_loading(reaction)
            assert_allclose(beam.f_loading[-1][-1], Rb1[i], rtol=0.01, err_msg='Reaction force does not match')

