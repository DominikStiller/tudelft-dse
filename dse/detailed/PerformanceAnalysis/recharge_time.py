import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from PIL import Image



if __name__ == '__main__':
    marsMap = np.asarray(Image.open("./Figures/Mars_map_BW.jpg"))
    plt.imshow(marsMap, interpolation='none', extent=[-180, 180, -90, 90])

    averageIrradiance = 590 * 132 * 2 / np.pi
    lat = np.radians(np.linspace(-90, 90, 180))
    lon = np.radians(np.linspace(-180, 180, 360))
    intensity = averageIrradiance * np.reshape(np.cos(lat), (lat.size, 1)) * np.ones(lon.size)

    recharge_time = 875 * 903 / (0.32 * intensity)
    for i in range(np.shape(recharge_time)[1]):
        recharge_time[0, i] = np.nan
        recharge_time[-1, i] = np.nan

    divNorm = mcolors.TwoSlopeNorm(vcenter=7.5, vmax=10, vmin=5)
    plt.imshow(recharge_time / 12, alpha=0.6, norm=divNorm, cmap='RdYlGn_r',
               interpolation='none', extent=[-180, 180, -90, 90])
    plt.colorbar(fraction=0.025, pad=0.04, label='Recharge time [sols]')
    # plt.hlines(np.arange(-90, 90, 18), xmin=-180, xmax=180, colors='black', alpha=0.25)
    plt.grid(axis='y', color='black', linestyle='--')
    plt.ylabel('Latitude [deg]')
    plt.xlabel('Longitude [deg]')
    plt.tight_layout()
    plt.savefig("./Figures/Recharge_time.png")
    plt.show()
