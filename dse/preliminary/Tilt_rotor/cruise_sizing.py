import numpy as np
import scipy.stats as sstats

from .constants import const


def update_wingspan(Print=True):
    global aircraftParameters
    aircraftParameters["wingspan"] = 3 * aircraftParameters["rotorRadius"]
    aircraftParameters["chord"] = aircraftParameters["wingArea"] / aircraftParameters["wingspan"]
    aircraftParameters["AR"] = aircraftParameters["wingspan"] / aircraftParameters["chord"]

    if aircraftParameters["AR"] < aircraftParameters["ARmin"]:
        # Try to decrease wingspan
        if Print:
            print(f"Chord is too large for the required AR")
        aircraftParameters["wingspan"] = np.sqrt(
            aircraftParameters["wingArea"] * aircraftParameters["ARmin"]
        )
        aircraftParameters["chord"] = (
            aircraftParameters["wingArea"] / aircraftParameters["wingspan"]
        )
        if aircraftParameters["wingspan"] < 2 * aircraftParameters["rotorRadius"]:
            # Didn't work, increase wingspan
            aircraftParameters["wingspan"] = 2 * aircraftParameters["rotorRadius"]
            aircraftParameters["chord"] = (
                aircraftParameters["wingArea"] / aircraftParameters["wingspan"]
            )
            aircraftParameters["AR"] = aircraftParameters["wingspan"] / aircraftParameters["chord"]
            while aircraftParameters["AR"] < aircraftParameters["ARmin"]:
                aircraftParameters["wingspan"] += 0.1
                aircraftParameters["chord"] = (
                    aircraftParameters["wingArea"] / aircraftParameters["wingspan"]
                )
                aircraftParameters["AR"] = (
                    aircraftParameters["wingspan"] / aircraftParameters["chord"]
                )
    if Print:
        print(f'Wing area = {aircraftParameters["wingArea"]}[m2]')
        print(f'Wingspan = {aircraftParameters["wingspan"]}[m]')
        print(f'Chord = {aircraftParameters["chord"]}[m]')
        print(f'AR = {aircraftParameters["wingspan"] / aircraftParameters["chord"]}')


def area(cl, MTOM, q, Print=True, changeConstants=None):
    if changeConstants is not None:
        const[changeConstants[0]] = changeConstants[1]

    g_mars = const["gravityMars"]

    aircraftParameters["wingArea"] = MTOM * g_mars / (cl * q)
    update_wingspan(Print)


def max_tipSpeed(cruiseVelocity):
    maxTipSpeed = np.sqrt((0.92 * 220) ** 2 - cruiseVelocity**2)
    tipSpeedVsCruiseSpeed = sstats.linregress(cruiseVelocity, maxTipSpeed)
    return tipSpeedVsCruiseSpeed.slope, tipSpeedVsCruiseSpeed.intercept
