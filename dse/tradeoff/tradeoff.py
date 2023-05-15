import numpy as np
from matplotlib import pyplot as plt

from dse.plotting import format_plot, save_plot
from dse.tradeoff.io import load_sheets


def calculate_scores(dfs, df_weights):
    expected_scores = []
    worst_scores = []
    best_scores = []

    for df in dfs:
        df["expected_weighted_score"] = df["expected_score"] * df.index.map(
            lambda c: df_weights.loc[c]["weight"]
        )
        df["worst_weighted_score"] = df["worst_score"] * df.index.map(
            lambda c: df_weights.loc[c]["weight"]
        )
        df["best_weighted_score"] = df["best_score"] * df.index.map(
            lambda c: df_weights.loc[c]["weight"]
        )

        total_expected = df["expected_weighted_score"].sum()
        total_worst = df["worst_weighted_score"].sum()
        total_best = df["best_weighted_score"].sum()

        expected_scores.append(total_expected)
        worst_scores.append(total_worst)
        best_scores.append(total_best)

    return expected_scores, worst_scores, best_scores


if __name__ == "__main__":
    design_names = ["Blended wing", "Conventional aircraft", "Tilt-rotor", "Multicopter", "Airship"]
    loaded_sheets = load_sheets("data/tradeoff.xlsx", design_names)
    expected_scores, worst_scores, best_scores = calculate_scores(*loaded_sheets)

    print("Scores:")
    for design_name, expected, worst, best in zip(
        design_names, expected_scores, worst_scores, best_scores
    ):
        print(f"  - {design_name}: {expected:.1f} ({worst:.1f} – {best:.1f})")

    print(f"Best design is {design_names[np.argmax(expected_scores)]}")

    fig, ax = plt.subplots(figsize=(7, 4))

    ax.scatter(design_names, expected_scores, marker="_", s=700, label="Expected")
    ax.scatter(design_names, worst_scores, marker="_", s=700, label="Worst case")
    ax.scatter(design_names, best_scores, marker="_", s=700, label="Best case")

    ax.legend()

    format_plot()
    save_plot(".", "tradeoff_results", type="png")
    plt.show()
