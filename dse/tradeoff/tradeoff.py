import numpy as np
from matplotlib import pyplot as plt

from dse.plotting import format_plot, save_plot
from dse.tradeoff.io import load_sheets


# df["expected_score"].apply(lambda s: df_scoring.loc[s]["score"])


def get_score(criterion, criterion_value, score_categories, df_scoring):
    score_categories = score_categories[criterion]

    first_closed_interval_criterion = list(score_categories.values())[1]
    ascending = first_closed_interval_criterion[1] > first_closed_interval_criterion[0]

    criterion_category = None
    for category, bounds in score_categories.items():
        if ascending:
            lower, upper = bounds
        else:
            upper, lower = bounds

        if lower is None:
            if criterion_value < upper:
                criterion_category = category
                break
        elif upper is None:
            if lower <= criterion_value:
                criterion_category = category
                break
        else:
            if lower <= criterion_value < upper:
                criterion_category = category
                break

    if not criterion_category:
        raise "Invalid criterion value"

    score = df_scoring.loc[criterion_category].iloc[0]
    return score


def calculate_total_score(
    df, df_weights, score_categories, df_scoring, maximum_score, column="expected"
):
    total_score = 0

    for criterion, criterion_value in df[column].items():
        criterion_score = get_score(criterion, criterion_value, score_categories, df_scoring)
        weight = df_weights.loc[criterion].iloc[0]
        total_score += weight * criterion_score

    return 100 * total_score / maximum_score


def calculate_score_regions(df_scoring, df_weights):
    regions = []

    total_weight = df_weights["weight"].sum()

    for score in df_scoring["score"].iloc[:-1].values:
        regions.append([None, total_weight * score])

    for i, upper in enumerate(regions[1:]):
        regions[i][0] = regions[i + 1][1]

    regions[-1][0] = 0

    return regions, regions[0][1]


if __name__ == "__main__":
    design_names = ["Blended wing", "Conventional aircraft", "Tilt-rotor", "Multicopter", "Airship"]
    dfs, df_weights, df_scoring, score_categories = load_sheets("data/tradeoff.xlsx", design_names)
    score_regions, maximum_score = calculate_score_regions(df_scoring, df_weights)
    total_scores = [
        calculate_total_score(df, df_weights, score_categories, df_scoring, maximum_score)
        for df in dfs
    ]

    # print("Scores:")
    # for design_name, expected, worst, best in zip(
    #     design_names, expected_scores, worst_scores, best_scores
    # ):
    #     print(f"  - {design_name}: {expected:.1f} ({worst:.1f} – {best:.1f})")
    #
    # print(f"Best design is {design_names[np.argmax(expected_scores)]}")

    fig, ax = plt.subplots(figsize=(7, 4))

    ax.scatter(design_names, total_scores, marker="_", s=700, label="Expected")
    # ax.scatter(design_names, worst_scores, marker="_", s=700, label="Worst case")
    # ax.scatter(design_names, best_scores, marker="_", s=700, label="Best case")

    for (lower, upper), color in zip(score_regions, ["375623", "70AD47", "FFC000", "ED7D31"]):
        lower = 100 * lower / maximum_score
        upper = 100 * upper / maximum_score
        ax.axhspan(lower, upper, color=f"#{color}", alpha=0.15)

    ax.set_ylim([0, 100])

    ax.legend()

    format_plot()
    save_plot(".", "tradeoff_results")
    plt.show()
