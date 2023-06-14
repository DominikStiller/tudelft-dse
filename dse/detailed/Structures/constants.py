import pandas as pd

const = {
    'rotorRadius': 8.2,
    'bladePerRotor': 8,
    'numRotors': 4,
    'MTOM': 3000,
    'g': 3.71,
    'rpm': 232 * 2 * 3.141592 / 60,
    'Airfoil': pd.read_csv("../../../tests/detailed/structures/S1223.dat", delimiter="\s+", dtype=float, skiprows=1, names=["x", "z"]),
    'cutout': 0.15,
    'fuselageHeight': 1.8,
    'engineMass': 13.2+21.9,
    'tailPoleLength': 15 - 2.1439,
    'tailWingLift': 0.1,
    'liftToDragTail': 15,
    'extraForceTail': 500,
    'q': 0.5 * 0.01 * 112 ** 2
}
