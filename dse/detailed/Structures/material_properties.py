import os

import pandas as pd


class Material:
    def __init__(
        self,
        density,
        E,
        G,
        yield_stress,
        thermal_exp_coefficient,
        tensile_strength,
        compressive_strength,
        shear_strength,
        poisson_ratio,
    ):
        self.rho = density
        self.E = E
        self.G = G
        self.sigmay = yield_stress
        self.a = thermal_exp_coefficient
        self.tensile = tensile_strength
        self.compressive = compressive_strength
        self.tau = shear_strength
        self.poisson = poisson_ratio


def read_material_properties(filename):
    material_properties = pd.read_csv(filename, delimiter=";", index_col=0)
    return material_properties


print(os.getcwd())
a = read_material_properties(os.path.normpath("dse/detailed/Structures/materials.csv"))

materials = dict()
for i, mat in enumerate(a.index):
    materials[mat.strip()] = Material(
        E=a["E"][i] * 1e9,
        G=a["G"][i] * 1e9,
        density=a["rho"][i],
        yield_stress=a["yield strength"][i] * 1e6,
        tensile_strength=a["tensile strength"][i] * 1e6,
        compressive_strength=a["compressive strength"][i] * 1e6,
        thermal_exp_coefficient=None,
        shear_strength=a["shear strength"][i] * 1e6,
        poisson_ratio=a["poisson"][i],
    )
