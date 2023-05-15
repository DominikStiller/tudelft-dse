import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from constants import const, aircraftParameters
from PIL import Image
import numpy as np



def v_i(thrust, radius):
    rho = const['airDensity']
    #N_rotor = const['RotorNumber']
    return np.sqrt(thrust/(rho*np.pi*radius**2))


def power(thrust, radius):
    return thrust * v_i(thrust, radius)


def size_power_subsystem(rotorRadius, takeOffThrust, cruiseThrust, cruiseTime, takeOffTime, surfaceArea, plot=False):
    # Define constants
    averageIrradiance = 2 / np.pi * const['irradiance']  # average intensity during the day
    cruiseTime /= 3600
    takeOffTime /= 3600
    collectingArea = surfaceArea / 2

    # Define the initial arrays
    lat = np.radians(np.linspace(-90, 90, 360))
    lon = np.radians(np.linspace(0, 360, 720))

    intensity = averageIrradiance * np.reshape(np.cos(lat), (lat.size, 1)) * np.ones(lon.size)

    # Power calculations
    takeOffPower = aircraftParameters['totalRotors'] * power(takeOffThrust/aircraftParameters['totalRotors'], rotorRadius)
    cruisePower = aircraftParameters['totalRotors'] * power(cruiseThrust/aircraftParameters['totalRotors'], rotorRadius)
    takeOffEnergy = takeOffPower * takeOffTime
    cruiseEnergy = cruisePower * cruiseTime

    powerSurplus = intensity * collectingArea - cruisePower
    rechargeCapability = np.copy(powerSurplus)
    rechargeCapability[np.where(rechargeCapability < 0)] = np.nan

    # Calculate in-flight recharge parameters
    rechargeTime = takeOffEnergy / rechargeCapability
    percentRechargeable = np.count_nonzero(~np.isnan(rechargeCapability)) / np.size(rechargeCapability)
    print(f'The minimum in-flight recharge time is {np.nanmin(rechargeTime)} hours')
    print(f'The average in-flight recharge time is {np.nanmean(rechargeTime)} hours')
    print(f'Recharge is possible in {100 * percentRechargeable}% of the planet\n')

    # Size the batteries and the arrays
    energyConsumption = 2 * takeOffEnergy + cruiseEnergy
    batteryMass = max(energyConsumption / const['batteryEnergyDensity'], takeOffPower/const['batteryPowerDensity'])
    cruiseBattery = cruiseEnergy / const['batteryEnergyDensity']
    panelMass = collectingArea * const['solarPanelDensity']
    print(f'Mass of the batteries = {batteryMass} kg')
    print(f'For cruise we need: {cruiseBattery*100/batteryMass}%')
    print(f'Mass of the solar panels = {panelMass} kg')

    if plot:
        # Plot the map of Mars
        marsMap = np.asarray(Image.open('../Figures/MarsMapBlackAndWhite.jpg'))
        plt.imshow(marsMap)

        # Plot the power surplus
        divNorm = mcolors.TwoSlopeNorm(vmin=np.min(powerSurplus)/1000, vcenter=0, vmax=np.max(powerSurplus)/1000)
        plt.pcolormesh(powerSurplus / 1000, alpha=0.6, cmap='RdYlGn', norm=divNorm)
        # plt.pcolormesh(powerSurplus / 1000, alpha=0.6, cmap='RdYlGn')
        plt.colorbar(fraction=0.025, pad=0.04, label='Energy surplus during cruise [kW]')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig('../Figures/EnergySurplusCruise.png')
        plt.show()

    return batteryMass, panelMass, powerSurplus


if __name__ == '__main__':
    # Define constants
    takeOffThrust = 1.1 * 3000 * 3.71
    cruiseThrust = 816.2
    rotorRadius = 14.1
    cruiseTime = 1e6/154
    takeOffTime = 600
    collectingArea = 60

    size_power_subsystem(rotorRadius, takeOffThrust, cruiseThrust, cruiseTime, takeOffTime, collectingArea, plot=True)
