# Sensitivities based on random disturbances
import matplotlib.pyplot as plt
import numpy as np

from dse.plotting import format_plot, save_plot
from dse.tradeoff.io import load_sheets
from dse.tradeoff.tradeoff import calculate_total_score, calculate_score_regions


def disturb_weights(dfs, df_weights, df_scoring, score_categories, maximum_score):
    N = 100
    std_relative = 0.1  # 10 % std

    scores_per_design = [[] for _ in design_names]

    for _ in range(N):
        df_weights_disturbed = df_weights.copy()
        dfs = [df.copy() for df in dfs]

        std = df_weights_disturbed["weight"] * std_relative
        df_weights_disturbed["weight"] += np.random.normal(scale=std)

        for i, df in enumerate(dfs):
            scores_per_design[i].append(
                calculate_total_score(
                    df, df_weights_disturbed, score_categories, df_scoring, maximum_score
                )
            )

    return scores_per_design


def disturb_scoring(dfs, df_weights, df_scoring, score_categories, maximum_score):
    N = 100
    std_relative = 0.1  # 10 % std

    scores_per_design = [[] for _ in design_names]

    for _ in range(N):
        df_scoring_disturbed = df_scoring.copy()
        dfs = [df.copy() for df in dfs]

        std = df_scoring_disturbed["score"] * std_relative
        df_scoring_disturbed["score"] += np.random.normal(scale=std)

        for i, df in enumerate(dfs):
            scores_per_design[i].append(
                calculate_total_score(
                    df, df_weights, score_categories, df_scoring_disturbed, maximum_score
                )
            )

    return scores_per_design


if __name__ == "__main__":
    design_names = ["Blended wing", "Conventional aircraft", "Tilt-rotor", "Multicopter", "Airship"]
    loaded_sheets = load_sheets("data/tradeoff.xlsx", design_names)
    dfs, df_weights, df_scoring, score_categories = loaded_sheets
    score_regions, maximum_score = calculate_score_regions(df_scoring, df_weights)

    scores_disturbed_scoring = disturb_scoring(*loaded_sheets, maximum_score)
    scores_disturbed_weights = disturb_weights(*loaded_sheets, maximum_score)

    fix, (ax_scoring, ax_weights) = plt.subplots(1, 2, figsize=(10, 5), sharey="all", sharex="all")

    ax_scoring.set_title("With scoring disturbed ($\\sigma$ = 10%)")
    ax_scoring.boxplot(scores_disturbed_scoring, labels=design_names)
    ax_scoring.set_ylabel("Total weighted score")

    ax_weights.set_title("With weights disturbed ($\\sigma$ = 10%)")
    ax_weights.boxplot(scores_disturbed_weights, labels=design_names)

    ax_scoring.set_xticklabels(ax_scoring.get_xticklabels(), rotation=15, ha="right")
    ax_weights.set_xticklabels(ax_weights.get_xticklabels(), rotation=15, ha="right")

    for (lower, upper), color in zip(score_regions, ["375623", "70AD47", "FFC000", "ED7D31"]):
        lower = 100 * lower / maximum_score
        upper = 100 * upper / maximum_score
        ax_scoring.axhspan(lower, upper, color=f"#{color}", alpha=0.15)
        ax_weights.axhspan(lower, upper, color=f"#{color}", alpha=0.15)

    ax_scoring.set_ylim([0, 100])

    format_plot()
    save_plot(".", "tradeoff_sensitivity_analysis")
    plt.show()
