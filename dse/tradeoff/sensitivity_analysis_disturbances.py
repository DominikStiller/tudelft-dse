# Sensitivities based on random disturbances
import matplotlib.pyplot as plt
import numpy as np

from dse.plotting import format_plot
from dse.tradeoff.io import load_sheets
from dse.tradeoff.tradeoff import calculate_scores


def disturb_weights(dfs, df_weights, df_scoring):
    N = 100
    std_relative = 0.1  # 10 % std

    scores_per_design = [[] for _ in design_names]

    for _ in range(N):
        df_weights_disturbed = df_weights.copy()
        dfs_disturbed = [df.copy() for df in dfs]

        std = df_weights["weight"] * std_relative
        df_weights_disturbed["weight"] += np.random.normal(scale=std)

        expected_scores, _, _ = calculate_scores(dfs_disturbed, df_weights_disturbed, df_scoring)
        for i, score in enumerate(expected_scores):
            scores_per_design[i].append(score)

    return scores_per_design


def disturb_values(dfs, df_weights, df_scoring):
    N = 100
    std_relative = 0.1  # 10 % std

    scores_per_design = [[] for _ in design_names]

    for _ in range(N):
        df_weights_disturbed = df_weights.copy()
        dfs_disturbed = []

        for df in dfs:
            std = np.abs(df["expected_score"].to_numpy()) * std_relative
            df_disturbed = df.copy()
            df_disturbed["expected_score"] += np.random.normal(scale=std)
            dfs_disturbed.append(df_disturbed)

        expected_scores, _, _ = calculate_scores(dfs_disturbed, df_weights_disturbed, df_scoring)
        for i, score in enumerate(expected_scores):
            scores_per_design[i].append(score)

    return scores_per_design


if __name__ == "__main__":
    # design_names = ["Blended wing", "Conventional aircraft", "Tilt-rotor", "Multicopter", "Airship"]
    design_names = ["Multicopter", "Airship"]
    loaded_sheets = load_sheets("data/tradeoff.xlsx", design_names)

    scores_expected, _, _ = calculate_scores(*loaded_sheets)
    scores_disturbed_values = disturb_values(*loaded_sheets)
    scores_disturbed_weights = disturb_weights(*loaded_sheets)

    fix, (ax_values, ax_weights) = plt.subplots(2, 1, figsize=(10, 10))

    ax_values.set_title("With scores for each criterion disturbed")
    ax_values.boxplot(scores_disturbed_values, labels=design_names)

    ax_weights.set_title("With weights disturbed")
    ax_weights.boxplot(scores_disturbed_weights, labels=design_names)

    format_plot()
    plt.show()
