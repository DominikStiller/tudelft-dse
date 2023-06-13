from numpy.testing import assert_allclose
from dse.detailed.Structures.main_structures import *
from dse.detailed.Structures.StructureClasses import Beam
from unittest import TestCase
import numpy as np


class TestBeam(TestCase):
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
        E = 337.5 * 1e9
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
        ...

    def test_size_wing(self):
        ...

    def test_size_tail(self):
        ...
