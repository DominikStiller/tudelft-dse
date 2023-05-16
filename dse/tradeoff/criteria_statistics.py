import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

from dse.plotting import save_plot
from dse.tradeoff.io import load_sheets
from dse.tradeoff.tradeoff import get_score


def plot_statistics(dfs):
    df = pd.DataFrame(
        {design_name: dfs[i]["expected"].to_numpy() for i, design_name in enumerate(design_names)},
        index=dfs[0].index,
    ).astype("float")

    criteria_names = list(dfs[0].index)

    fig, ax = plt.subplots(figsize=(10, 5), tight_layout=True)
    im = ax.imshow(df.to_numpy().T / df.to_numpy().mean(axis=1)[:, None].T, cmap="Reds")
    fig.colorbar(im)

    std_pct = df.std(axis=1) / df.mean(axis=1) * 100

    _add_labels(ax, criteria_names, std_pct)

    save_plot(".", "tradeoff_values", type="png")
    plt.show()

    with pd.option_context(
        "display.max_rows", None, "display.max_columns", None, "display.width", None
    ):
        print(df)


def plot_scores(dfs, df_scoring, score_categories):
    df = pd.DataFrame(
        {design_name: dfs[i]["expected"].to_numpy() for i, design_name in enumerate(design_names)},
        index=dfs[0].index,
    ).astype("float")

    criteria_names = list(dfs[0].index)

    df_scores = df.copy()
    for design_name in design_names:
        for criterion, criterion_value in df_scores[design_name].items():
            df_scores[design_name].loc[criterion] = get_score(
                criterion, criterion_value, score_categories, df_scoring
            )

    fig, ax = plt.subplots(figsize=(10, 5), tight_layout=True)
    im = ax.imshow(
        df_scores.values.T,
        cmap=ListedColormap(["#C00000", "#ED7D31", "#FFC000", "#70AD47", "#375623"]),
    )
    fig.colorbar(im)

    _add_labels(ax, criteria_names)

    plt.show()


def _add_labels(ax, criteria_names, std_pct=None):
    ax.set_yticks(np.arange(len(design_names)), labels=design_names)
    if std_pct is None:
        ax.set_xticks(
            np.arange(len(criteria_names)),
            labels=[f"{criterion_name}" for j, criterion_name in enumerate(criteria_names)],
        )
    else:
        ax.set_xticks(
            np.arange(len(criteria_names)),
            labels=[
                f"{criterion_name} ($\sigma$ = {std_pct.iloc[j]:.0f} %)"
                for j, criterion_name in enumerate(criteria_names)
            ],
        )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=20, ha="right")

    for i in range(len(design_names)):
        for j, criterion_name in enumerate(criteria_names):
            ax.text(
                j,
                i,
                f"{dfs[i].loc[criterion_name]['expected']:.1f}",
                ha="center",
                va="center",
                weight="bold" if criterion_name in selected_criteria else "normal",
            )


if __name__ == "__main__":
    design_names = ["Blended wing", "Conventional aircraft", "Tilt-rotor", "Multicopter", "Airship"]

    dfs, df_weights, _, _ = load_sheets("data/tradeoff.xlsx", design_names, selected_only=False)
    selected_criteria = df_weights[df_weights["selected"] == "x"].index

    plot_statistics(dfs)

    dfs, df_weights, df_scoring, score_categories = load_sheets(
        "data/tradeoff.xlsx", design_names, selected_only=True
    )
    plot_scores(dfs, df_scoring, score_categories)
