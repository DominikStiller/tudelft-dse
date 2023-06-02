import matplotlib.pyplot as plt

from dse.plotting import format_plot, save_plot

if __name__ == "__main__":
    # Expected
    kg_per_h = 130.16
    takeoff_plus_land = 130.16 / 2

    # # Best
    # kg_per_h = 44.1
    # takeoff_plus_land = 67.8

    # # Worst
    # kg_per_h = 145.78
    # takeoff_plus_land = 145.78 / 2

    km_per_kg = (400 / kg_per_h) / 1.3
    mtom = 2700  # kg
    oem = 1747  # kg

    def payload2fuel(payload):
        return mtom - oem - payload

    def fuel2payload(fuel):
        return mtom - oem - fuel

    max_payload = mtom - oem - takeoff_plus_land
    print(f"{max_payload=}")

    max_payload_at_1000km = mtom - oem - takeoff_plus_land - 1000 / km_per_kg
    print(f"{max_payload_at_1000km=}")

    fuel_mass_at_350kg_payload = mtom - oem - 350
    max_range_at_350kg_payload = km_per_kg * (fuel_mass_at_350kg_payload - takeoff_plus_land)
    print(f"{fuel_mass_at_350kg_payload=}")
    print(f"{max_range_at_350kg_payload=}")

    max_fuel_at_max_payload = mtom - oem - takeoff_plus_land - max_payload
    max_range_at_max_payload = max_fuel_at_max_payload * km_per_kg
    print(f"{max_range_at_max_payload=}")

    max_fuel_at_no_payload = mtom - oem - takeoff_plus_land
    max_range_at_no_payload = max_fuel_at_no_payload * km_per_kg
    print(f"{max_range_at_no_payload=}")

    fig, ax = plt.subplots(figsize=(9, 2.5))
    secax = ax.secondary_yaxis("right", functions=(payload2fuel, fuel2payload))
    secax.set_ylabel("Fuel mass [kg]")

    ax.plot(
        [0, max_range_at_max_payload, max_range_at_no_payload],
        [max_payload, max_payload, 0],
        c="black",
    )
    ax.vlines(1000, 0, max_payload_at_1000km, color="black", ls=":")
    ax.hlines(350, 0, max_range_at_350kg_payload, color="black", ls=":")
    ax.set_xlabel("Range [km]")
    ax.set_ylabel("Payload mass [kg]")
    ax.set_xlim([0, max_range_at_no_payload * 1.05])
    ax.set_ylim([0, max_payload * 1.05])

    format_plot()
    save_plot(".", "multicopter_payload_range")
    plt.show()
