import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from constants import const, aircraftParameters
from PIL import Image
import numpy as np



def v_i(thrust, radius):
    rho = const['airDensity']
    return np.sqrt(thrust/(2*rho*np.pi*radius**2))


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
    takeOffPower = aircraftParameters['totalRotors'] * power(takeOffThrust, rotorRadius)
    # cruisePower = aircraftParameters['totalRotors'] * power(cruiseThrust/aircraftParameters['totalRotors'], rotorRadius)
    cruisePower = cruiseThrust * const['cruiseSpeed']
    takeOffEnergy = takeOffPower * takeOffTime
    cruiseEnergy = cruisePower * cruiseTime

    powerSurplus = intensity * collectingArea - cruisePower
    rechargeCapability = np.copy(powerSurplus)
    rechargeCapability[np.where(rechargeCapability < 0)] = np.nan

    # Calculate in-flight recharge parameters
    try:
        percentRechargeable = np.count_nonzero(~np.isnan(rechargeCapability)) / np.size(rechargeCapability)
    except:
        rechargeTime = (2 * takeOffEnergy + cruiseEnergy) / (intensity * collectingArea)
        print(f'We cannot recharge in flight - average recharge time landed is {np.mean(rechargeTime)}')

    rechargeTime = (2 * takeOffEnergy + cruiseEnergy) / (intensity * collectingArea) / 3600
    print(f'Cruise power = {cruisePower} W')
    print(f'Cruise energy = {cruiseEnergy} Wh')
    print(f'Take off power = {takeOffPower} W')
    print(f'Take off energy = {takeOffEnergy} Wh')
    print(f'The minimum in-flight recharge time is {np.nanmin(rechargeTime)} hours')
    print(f'The average in-flight recharge time is {np.nanmean(rechargeTime)} hours')
    print(f'Average recharge time landed is {np.mean(rechargeTime[1, -2])}')
    print(f'Recharge is possible in {100 * percentRechargeable}% of the planet while in cruise\n')

    # Size the batteries and the arrays
    energyConsumption = 2 * takeOffEnergy + cruiseEnergy
    batteryMass = max(2*takeOffEnergy / const['takeoffBatteryEnergyDensity'], takeOffPower/const['takeoffBatteryPowerDensity'])
    # cruiseBattery = cruiseEnergy / const['batteryEnergyDensity']
    cruiseBattery = cruiseEnergy / const['takeoffBatteryEnergyDensity']
    batteryMass += cruiseBattery
    panelMass = collectingArea * const['solarPanelDensity']

    # Apply safety margins
    batteryMass *= 1.3/0.95

    print(f'Mass of the batteries = {batteryMass} kg')
    print(f'Volume of the batteries = {energyConsumption/const["batteryVolume"]}')
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
        plt.colorbar(fraction=0.025, pad=0.04, label='Power surplus during cruise [kW]')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig('../Figures/PowerSurplusCruise.png')
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
