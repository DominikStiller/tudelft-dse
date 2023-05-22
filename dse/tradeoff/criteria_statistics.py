import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

from dse.plotting import save_plot
from dse.tradeoff.io import load_sheets
from dse.tradeoff.tradeoff import get_score


def plot_statistics(dfs, df_weights):
    df = pd.DataFrame(
        {design_name: dfs[i]["expected"].to_numpy() for i, design_name in enumerate(design_names)},
        index=dfs[0].index,
    ).astype("float")

    df = df.dropna()
    df_weights = df_weights.loc[df.index]

    fig, ax = plt.subplots(figsize=(10, 5), tight_layout=True)
    im = ax.imshow(
        (df.to_numpy() - np.nanmin(df.to_numpy(), axis=1)[:, None]).T
        / (np.nanmax(df.to_numpy(), axis=1) - np.nanmin(df.to_numpy(), axis=1))[:, None].T,
        cmap="Reds",
    )
    fig.colorbar(im)

    std_pct = df.std(axis=1) / df.mean(axis=1) * 100

    _add_labels(ax, df_weights, std_pct=std_pct)

    fig.tight_layout(pad=0.1, h_pad=0.4, w_pad=0.4)
    save_plot("tradeoff", "tradeoff_values")
    plt.show()

    with pd.option_context(
        "display.max_rows", None, "display.max_columns", None, "display.width", None
    ):
        print(df)


def plot_summary(dfs, df_weights, df_scoring, score_categories):
    df = pd.DataFrame(
        {design_name: dfs[i]["expected"].to_numpy() for i, design_name in enumerate(design_names)},
        index=dfs[0].index,
    )

    df_scores = df.copy()
    for design_name in design_names:
        for criterion, criterion_value in df_scores[design_name].items():
            df_scores[design_name].loc[criterion] = get_score(
                criterion, criterion_value, score_categories, df_scoring
            )

    scores = df_scores.values.T

    fig, ax = plt.subplots(figsize=(10, 5), tight_layout=True)
    im = ax.imshow(
        scores,
        cmap=ListedColormap(["#C00000", "#ED7D31", "#FFC000", "#70AD47", "#548235"]),
    )
    cbar = fig.colorbar(im, boundaries=[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], ticks=[0, 1, 2, 3, 4])
    cbar.ax.set_yticklabels(df_scoring.index[::-1])

    _add_labels(
        ax,
        df_weights,
        selected_bold=False,
    )

    fig.tight_layout(pad=0.1, h_pad=0.4, w_pad=0.4)
    save_plot("tradeoff", "tradeoff_summary")
    plt.show()


def _add_labels(ax, df_weights, selected_bold=True, std_pct=None):
    criteria_names = list(df_weights.index)
    units = df_weights["unit"].fillna("-")

    ax.set_yticks(np.arange(len(design_names)), labels=design_names)
    if std_pct is None:
        ax.set_xticks(
            np.arange(len(criteria_names)),
            labels=[
                f"{criterion_name} [{unit}]"
                for j, (criterion_name, unit) in enumerate(zip(criteria_names, units))
            ],
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
                weight="bold"
                if (selected_bold and criterion_name in selected_criteria)
                else "normal",
            )


if __name__ == "__main__":
    design_names = ["Blended wing", "Conventional aircraft", "Tilt-rotor", "Multicopter", "Airship"]

    dfs, df_weights, _, _ = load_sheets("data/tradeoff.xlsx", design_names, selected_only=False)
    selected_criteria = df_weights[df_weights["selected"] == "x"].index

    plot_statistics(dfs, df_weights)

    dfs, df_weights, df_scoring, score_categories = load_sheets(
        "data/tradeoff.xlsx", design_names, selected_only=True
    )
    plot_summary(dfs, df_weights, df_scoring, score_categories)
