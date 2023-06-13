from numpy.testing import assert_allclose
from dse.detailed.Structures.main_structures import *
from dse.detailed.Structures.StructureClasses import Beam
from unittest import TestCase
import numpy as np


class TestMainStructures(TestCase):
    def test_rotor_vibrations(self):
        # According to literature:
        # Fixed-free: wn = 1/2pi * (3.5156/L^2) * sqrt(EI/rho_L)

        # Square beam with 5mm thickness
        # Analytical solution
        x = 0.1
        z = 0.1
        t = 0.005
        A = x*z - (x-t)*(z-t)

        L = 2
        rho_L = A * 1645
        E = 337.5 * 1e9
        I = x*t**3/6 + 2*x*t*(z-t)**2/4 + t*x**3/6
        wn0 = (3.5156 / L**2) * np.sqrt(E * I / rho_L)

        # Code solution
        beam = Beam(
            width=x,
            height=z,
            length=L,
            cross_section='square',
            material='CFRCy',
            fixing_points=np.array([[x/2], [z/2]])*np.ones(100)
        )

        beam.Bi = A / np.shape(beam.x)[0] * np.ones(np.shape(beam.x))
        wn = rotor_vibrations(beam, reinforce=False, overwrite_I=I)[0][0]

        assert_allclose(wn, wn0, rtol=1e-2, err_msg='Natural frequency is incorrect')

    def test_wing_vibrations(self):
        # According to literature:
        # Fixed-pinned: wn = (15.4128/L^2) * sqrt(EI/rho_L)

        # Square beam with 5mm thickness
        # Analytical solution
        x = 0.1
        z = 0.1
        t = 0.005
        A = x * z - (x - t) * (z - t)

        L = 2
        rho_L = A * 1645
        E = 147.5 * 1e9
        I = x * t ** 3 / 6 + 2 * x * t * (z - t) ** 2 / 4 + t * x ** 3 / 6
        wn0 = (15.4128 / L ** 2) * np.sqrt(E * I / rho_L)

        # Code solution
        beam = Beam(
            width=x,
            height=z,
            length=L,
            cross_section='square',
            material='CFRCy',
            fixing_points=np.array([[x / 2], [z / 2]]) * np.ones(100)
        )

        parameters = np.array([E, I, 1645, A, L])
        wn = wing_vibrations(parameters)[0][0]

        assert_allclose(wn, wn0, rtol=1e-2, err_msg='Natural frequency is incorrect')

    def test_size_rotor_blades(self):
        s_max = 939 * 1e6 / 3
        front, rear, mr = size_rotor_blades(overwrite=True)
        self.assertTrue(np.all(np.abs(front.sigma) < s_max), msg='Stress exceeds allowable compressive stress')
        self.assertTrue(np.all(np.abs(rear.sigma) < s_max), msg='Stress exceeds allowable compressive stress')

    def test_size_wing(self):
        s_max = 1072.5 * 1e6 / 3
        span = np.flip(np.array([45, 40, 35, 30]))
        rootChord = np.flip(np.array([4, 4.3333, 4.9481, 5.78]))

        for i in range(len(span)):
            wing = size_wing(span[i], rootChord[i], 0.5, i)[0]
            self.assertTrue(np.all(np.abs(wing.sigma) < s_max), msg='Stress exceeds allowable compressive stress')

    def test_size_tail(self):
        s_max0 = 1072.5 * 1e6 / 3
        s_max1 = 185 * 1e6 / 1.5
        hs, vs, tp = size_tail()
        self.assertTrue(np.all(np.abs(hs.sigma) < s_max0), msg='Stress exceeds allowable compressive stress')
        self.assertTrue(np.all(np.abs(vs.sigma) < s_max1), msg='Stress exceeds allowable compressive stress')
