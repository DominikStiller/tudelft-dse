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


def calculate_total_score(df, df_weights, score_categories, df_scoring, column="expected"):
    total_score = 0

    for criterion, criterion_value in df[column].items():
        criterion_score = get_score(criterion, criterion_value, score_categories, df_scoring)
        weight = df_weights.loc[criterion].iloc[0]
        total_score += weight * criterion_score

    return total_score


if __name__ == "__main__":
    design_names = ["Blended wing", "Conventional aircraft", "Tilt-rotor", "Multicopter", "Airship"]
    dfs, df_weights, df_scoring, score_categories = load_sheets("data/tradeoff.xlsx", design_names)
    total_scores = [
        calculate_total_score(df, df_weights, score_categories, df_scoring) for df in dfs
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

    ax.legend()

    format_plot()
    save_plot(".", "tradeoff_results")
    plt.show()
