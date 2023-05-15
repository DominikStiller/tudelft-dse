# Sensitivities based on the ranges of uncertainty
from dse.tradeoff.io import load_sheets

if __name__ == "__main__":
    # design_names = ["Blended wing", "Conventional aircraft", "Tilt-rotor", "Multicopter", "Airship"]
    design_names = ["Multicopter", "Airship"]
    loaded_sheets = load_sheets("data/tradeoff.xlsx", design_names)

    scores_values = disturb_values(*loaded_sheets)

    fix, (ax_values, ax_weights) = plt.subplots(2, 1, figsize=(10, 10))

    ax_values.boxplot(scores_values, labels=design_names)
    ax_weights.boxplot(scores_weights, labels=design_names)

    format_plot()
    plt.show()
