from dse.preliminary.Tilt_rotor.RotorEngineSizing import RadiusMassElementMomentum
from dse.preliminary.Tilt_rotor.constants import *
from dse.preliminary.Tilt_rotor.cruise_sizing import *
from dse.preliminary.Tilt_rotor.power_sizing import size_power_subsystem


def RotorSizingValidation():
    referencedataV22 = {
        "rotorRadius": 11.6 / 2,
        "maxGrossWeight": 23859,
        "propulsionPower": 2749000 * 2,
        "NumberRotor": 2,
        "NumberBlade": 3,
        "TipSpeed": 0.9,
    }
    referencedataV280 = {
        "rotorRadius": 10.7 / 2,
        "maxGrossWeight": 17200,
        "propulsionPower": 4470000,
        "NumberRotor": 2,
        "NumberBlade": 3,
        "TipSpeed": 0.75,
    }
    referencedataAW609 = {
        "rotorRadius": 7.92 / 2,
        "maxGrossWeight": 7620,
        "propulsionPower": 1940 * 745.7,
        "NumberRotor": 2,
        "NumberBlade": 3,
        "TipSpeed": 0.7,
    }
    datas = [referencedataV22, referencedataV280, referencedataAW609]

    calcR, calcT, calcHP, calcPower, Rotor_mass = ([],) * 5
    (
        ReferenceRadius,
        ReferenceWeight,
        ReferencePower,
        ReferenceNrotor,
        ReferenceNblade,
        ReferenceTipSpeed,
    ) = ([],) * 6
    i = 0
    for elem in datas:
        ReferenceRadius = np.append(ReferenceRadius, elem["rotorRadius"])
        ReferenceWeight = np.append(ReferenceWeight, elem["maxGrossWeight"])
        ReferencePower = np.append(ReferencePower, elem["propulsionPower"])
        ReferenceNrotor = np.append(ReferenceNrotor, elem["NumberRotor"])
        ReferenceNblade = np.append(ReferenceNblade, elem["NumberBlade"])
        ReferenceTipSpeed = np.append(ReferenceTipSpeed, elem["TipSpeed"])

    for i in range(3):
        # Define constants for Earth:
        const["gravityMars"] = 9.80665
        const["airDensity"] = 1.225
        const["soundSpeed"] = 343
        R, T, hp, P, m = RadiusMassElementMomentum(
            M=ReferenceWeight[i],
            N_rotors=ReferenceNrotor[i],
            N_blades=ReferenceNblade[i],
            coaxial=False,
            V_tip=ReferenceTipSpeed[i] * const["soundSpeed"],
            print_results=False,
        )

        calcR = np.append(calcR, R)
        calcT = np.append(calcT, T)
        calcPower = np.append(calcPower, P)

    DeltaR = (ReferenceRadius - calcR) / ReferenceRadius
    DeltaPower = (ReferencePower - calcPower) / ReferencePower
    print(f"Difference in Rotor Radius of {DeltaR*100}[m]")
    print(f"Difference in power of {DeltaPower*100}[W]")


def AreaValidation():
    const["gravityMars"] = 9.80665
    aircraftParameters["rotorRadius"] = 1
    C_l = 1.2
    m = 100  # kg
    rho = 1.225  # kg/m3
    V = 100  # m/s
    p_dyn = 0.5 * rho * V * V
    area(C_l, m, p_dyn, Print=False, changeConstants=None)
    a = aircraftParameters["wingArea"]
    S = (m * 9.80665) / (0.5 * C_l * rho * V * V)
    print(f"Area function is valid: {S- a <= 1e-6}")


def DragValidation():
    print(
        "The drag estimations stem from Torenbeek's book: Synthesis of Subsonic Airplane Design.\n As such, the methods are deemed valid"
    )


def PowerSizingValidation():
    rotorRadius = 5
    rho = 0.01
    A = np.pi * rotorRadius**2
    takeOffThrust = 1000  # N
    cruiseThrust = 100  # N
    cruiseTime = 1.5  # hrs
    takeOffTime = 0.1  # hrs
    surfaceArea = 10  # m2
    const["cruiseSpeed"] = 100
    aircraftParameters["totalRotors"] = 1

    # Thrust Power:
    P_TO = (
        aircraftParameters["totalRotors"] * takeOffThrust * np.sqrt(takeOffThrust / (2 * rho * A))
    )
    E_TO = P_TO * takeOffTime

    # Cruise
    P_cruise = cruiseThrust * const["cruiseSpeed"]
    E_cruise = P_cruise * cruiseTime

    PowerDensity, EnergyDensity = (1317, 437)

    TO_battery = max(P_TO / PowerDensity, E_TO / EnergyDensity)
    Cruise_battery = 0
    if TO_battery * EnergyDensity >= E_cruise + 2 * E_TO:
        Cruise_battery = max(
            (E_cruise - (TO_battery * EnergyDensity - E_TO)) / EnergyDensity,
            P_cruise / PowerDensity,
        )

    RefBattery = (TO_battery + Cruise_battery) * 1.3 / 0.95

    calcBatteryMass, panelMass, powerSurplus = size_power_subsystem(
        rotorRadius,
        takeOffThrust,
        cruiseThrust,
        cruiseTime,
        takeOffTime,
        surfaceArea,
        plot=False,
        Print=False,
    )
    print(
        f"Difference between calculated and required battery mass: {100*(RefBattery-calcBatteryMass)/RefBattery}%"
    )


def PayloadRangeValidation():
    R = 1
    mass_rotors = 1
    Mass_design = 1
    N_ult = 1
    AR = 1
    wingbraced = 1
    V_cr = 1
    E_density = 1
    P_density_TO = 1
    E_density_TO = 1
    Mass_solar = 1
    minpayload = 1
    minRange = 1
    availableSpace = 1


RotorSizingValidation()
# AreaValidation()
# PowerSizingValidation()
