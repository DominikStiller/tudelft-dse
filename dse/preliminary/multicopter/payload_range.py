import matplotlib.pyplot as plt

from dse.plotting import format_plot, save_plot

if __name__ == "__main__":
    km_per_kg = (400 / 146) / 1.3
    mtom = 2700  # kg
    oem = 1747  # kg
    max_payload = 350  # kg
    takeoff_plus_land = 36.4 * 2  # kg

    max_fuel_at_max_payload = mtom - oem - takeoff_plus_land - max_payload
    max_range_at_max_payload = max_fuel_at_max_payload * km_per_kg
    print(max_range_at_max_payload)

    max_fuel_at_no_payload = mtom - oem
    max_range_at_no_payload = max_fuel_at_no_payload * km_per_kg
    print(max_range_at_no_payload)

    # # For flying wing
    # max_payload = 350
    # max_range_at_max_payload = 1000
    # max_range_at_no_payload = 6241.043

    fig, ax = plt.subplots(figsize=(10, 3))

    ax.plot(
        [0, max_range_at_max_payload, max_range_at_no_payload],
        [max_payload, max_payload, 0],
        c="black",
    )
    ax.vlines(max_range_at_max_payload, 0, max_payload, color="black", ls=":")
    ax.set_xlabel("Range [km]")
    ax.set_ylabel("Payload mass [kg]")
    ax.set_xlim([0, max_range_at_no_payload * 1.05])
    ax.set_ylim([0, max_payload * 1.05])

    format_plot()
    save_plot(".", "multicopter_payload_range")
    plt.show()
