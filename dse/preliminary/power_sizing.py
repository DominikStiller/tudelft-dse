import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np


def v_i(thrust, radius):
    rho = 0.01
    return np.sqrt(thrust/(2*rho*np.pi*radius**2))


def power(thrust, radius):
    return thrust * v_i(thrust, radius)


if __name__ == '__main__':
    # Define constants
    averageIrradiance = 2/np.pi * 590  # average intensity during the day
    takeOffThrust = 1.1 * 3000 * 3.71
    cruiseThrust = 816.2
    rotorRadius = 14.1
    cruiseTime = 1e6/154/3600
    takeOffTime = 1/6
    collectingArea = 60
    batteryEnergyDensity = 350
    solarPanelDensity = 1.76  # https://www.spectrolab.com/DataSheets/Panel/panels.pdf


    # Define the initial arrays
    lat = np.radians(np.arange(-90, 90, 0.5))
    lon = np.radians(np.arange(0, 360, 0.5))

    intensity = averageIrradiance * np.reshape(np.cos(lat), (lat.size, 1)) * np.ones(lon.size)

    # Power calculations
    takeOffPower = power(takeOffThrust, rotorRadius)
    cruisePower = power(cruiseThrust, rotorRadius)
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
    print(f'Recharge is possible in {100*percentRechargeable}% of the planet\n')

    # Size the batteries and the arrays
    energyConsumption = 2 * takeOffEnergy + cruiseEnergy
    batteryMass = energyConsumption / batteryEnergyDensity
    panelMass = collectingArea * solarPanelDensity
    print(f'Mass of the batteries = {batteryMass} kg')
    print(f'Mass of the solar panels = {panelMass} kg')


    # # Plot the map of Mars
    # marsMap = np.asarray(Image.open('Figures/MarsMapBlackAndWhite.jpg'))
    # plt.imshow(marsMap)
    #
    # # Plot the power surplus
    # divNorm = mcolors.TwoSlopeNorm(vmin=np.min(powerSurplus)/1000, vcenter=0, vmax=np.max(powerSurplus)/1000)
    # plt.pcolormesh(powerSurplus / 1000, alpha=0.6, cmap='RdYlGn', norm=divNorm)
    # plt.colorbar(fraction=0.025, pad=0.04, label='Energy surplus during cruise [kW]')
    # plt.axis('off')
    # plt.tight_layout()
    # plt.savefig('Figures/EnergySurplusCruise.png')
    # plt.show()
