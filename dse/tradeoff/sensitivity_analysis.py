# Sensitivities based on random perturbances
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullLocator

from dse.plotting import format_plot, save_plot
from dse.tradeoff.io import load_sheets
from dse.tradeoff.tradeoff import calculate_total_score, calculate_score_regions


def perturb_weights(dfs, df_weights, df_scoring, score_categories):
    N = 1000
    max_spread = 5

    scores_per_design = [[] for _ in design_names]

    for _ in range(N):
        df_weights_perturbed = df_weights.copy()
        dfs = [df.copy() for df in dfs]

        df_weights_perturbed["weight"] = np.random.uniform(1, max_spread, len(df_weights.index))

        _, maximum_score = calculate_score_regions(df_scoring, df_weights_perturbed)

        for i, df in enumerate(dfs):
            scores_per_design[i].append(
                calculate_total_score(
                    df, df_weights_perturbed, score_categories, df_scoring, maximum_score
                )
            )

    return scores_per_design


if __name__ == "__main__":
    np.random.seed(42)

    design_names = ["Blended wing", "Biplane aircraft", "Tiltrotor", "Multicopter", "Airship"]
    loaded_sheets = load_sheets("data/tradeoff.xlsx", design_names)
    dfs, df_weights, df_scoring, score_categories = loaded_sheets
    score_regions, maximum_score = calculate_score_regions(df_scoring, df_weights)

    scores_perturbed = perturb_weights(*loaded_sheets)

    print("Perturbed weights:")
    for design_name, scores in zip(design_names, scores_perturbed):
        median = np.median(scores)
        q1 = np.percentile(scores, 25)
        q3 = np.percentile(scores, 75)
        iqr = q3 - q1

        print(f"  - {design_name}: {median:.0f} ({q1:.0f} – {q3:.0f}, IQR = {iqr:.0f})")

    fix, ax = plt.subplots(figsize=(10, 3), sharey="all", sharex="all")

    ax.boxplot(scores_perturbed, labels=design_names)

    for (lower, upper), color in zip(score_regions, ["375623", "70AD47", "FFC000", "ED7D31"]):
        lower = 100 * lower / maximum_score
        upper = 100 * upper / maximum_score
        ax.axhspan(lower, upper, color=f"#{color}", alpha=0.3)

    ax.set_ylim([0, 100])
    ax.set_ylabel("Total weighted score")
    ax.set_yticks([0, 25, 50, 75, 100])
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=15, ha="right")

    format_plot(xlocator=NullLocator())
    save_plot("tradeoff", "tradeoff_sensitivity_analysis")
    plt.show()
