# Sensitivities based on random perturbances
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullLocator

from dse.plotting import format_plot, save_plot
from dse.tradeoff.io import load_sheets
from dse.tradeoff.tradeoff import calculate_total_score, calculate_score_regions


def perturb_weights(dfs, df_weights, df_scoring, score_categories, maximum_score):
    N = 1000
    std_relative = 0.1  # 10 % std

    scores_per_design = [[] for _ in design_names]

    for _ in range(N):
        df_weights_perturbed = df_weights.copy()
        dfs = [df.copy() for df in dfs]

        std = df_weights_perturbed["weight"] * std_relative
        df_weights_perturbed["weight"] += np.random.normal(scale=std)

        for i, df in enumerate(dfs):
            scores_per_design[i].append(
                calculate_total_score(
                    df, df_weights_perturbed, score_categories, df_scoring, maximum_score
                )
            )

    return scores_per_design


def perturb_scoring(dfs, df_weights, df_scoring, score_categories, maximum_score):
    N = 1000
    std_relative = 0.1  # 10 % std

    scores_per_design = [[] for _ in design_names]

    for _ in range(N):
        df_scoring_perturbed = df_scoring.copy()
        dfs = [df.copy() for df in dfs]

        std = df_scoring_perturbed["score"] * std_relative
        df_scoring_perturbed["score"] += np.random.normal(scale=std)

        for i, df in enumerate(dfs):
            scores_per_design[i].append(
                calculate_total_score(
                    df, df_weights, score_categories, df_scoring_perturbed, maximum_score
                )
            )

    return scores_per_design


if __name__ == "__main__":
    np.random.seed(42)

    design_names = ["Blended wing", "Conventional aircraft", "Tilt-rotor", "Multicopter", "Airship"]
    loaded_sheets = load_sheets("data/tradeoff.xlsx", design_names)
    dfs, df_weights, df_scoring, score_categories = loaded_sheets
    score_regions, maximum_score = calculate_score_regions(df_scoring, df_weights)

    scores_perturbed_scoring = perturb_scoring(*loaded_sheets, maximum_score)
    scores_perturbed_weights = perturb_weights(*loaded_sheets, maximum_score)

    print("Perturbed scoring:")
    for design_name, scores in zip(design_names, scores_perturbed_scoring):
        median = np.median(scores)
        q1 = np.percentile(scores, 25)
        q3 = np.percentile(scores, 75)
        iqr = q3 - q1

        print(f"  - {design_name}: {median:.0f} ({q1:.0f} – {q3:.0f}, IQR = {iqr:.0f})")

    print("Perturbed weights:")
    for design_name, scores in zip(design_names, scores_perturbed_weights):
        median = np.median(scores)
        q1 = np.percentile(scores, 25)
        q3 = np.percentile(scores, 75)
        iqr = q3 - q1

        print(f"  - {design_name}: {median:.0f} ({q1:.0f} – {q3:.0f}, IQR = {iqr:.0f})")

    fix, (ax_scoring, ax_weights) = plt.subplots(1, 2, figsize=(10, 5), sharey="all", sharex="all")

    ax_scoring.set_title("With scoring perturbed")
    ax_scoring.boxplot(scores_perturbed_scoring, labels=design_names)
    ax_scoring.set_ylabel("Total weighted score")

    ax_weights.set_title("With weights perturbed")
    ax_weights.boxplot(scores_perturbed_weights, labels=design_names)

    ax_scoring.set_xticklabels(ax_scoring.get_xticklabels(), rotation=15, ha="right")
    ax_weights.set_xticklabels(ax_weights.get_xticklabels(), rotation=15, ha="right")

    for (lower, upper), color in zip(score_regions, ["375623", "70AD47", "FFC000", "ED7D31"]):
        lower = 100 * lower / maximum_score
        upper = 100 * upper / maximum_score
        ax_scoring.axhspan(lower, upper, color=f"#{color}", alpha=0.15)
        ax_weights.axhspan(lower, upper, color=f"#{color}", alpha=0.15)

    ax_scoring.set_ylim([0, 100])
    ax_scoring.set_yticks([0, 25, 50, 75, 100])

    format_plot(xlocator=NullLocator())
    save_plot("tradeoff", "tradeoff_sensitivity_analysis")
    plt.show()
