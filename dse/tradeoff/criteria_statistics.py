import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from dse.plotting import save_plot
from dse.tradeoff.io import load_sheets


if __name__ == "__main__":
    design_names = ["Blended wing", "Conventional aircraft", "Tilt-rotor", "Multicopter", "Airship"]
    dfs, df_weights = load_sheets(
        "data/tradeoff.xlsx", design_names, selected_only=False, convert_score=False
    )
    selected_criteria = df_weights[df_weights["selected"] == "x"].index
    criteria_names = list(dfs[0].index)

    df = pd.DataFrame(
        {design_name: dfs[i]["expected"].to_numpy() for i, design_name in enumerate(design_names)},
        index=dfs[0].index,
    ).astype("float")
    # Remove qualitative criteria
    df = df.iloc[:-3]
    criteria_names = criteria_names[:-3]

    fig, ax = plt.subplots(figsize=(10, 5), tight_layout=True)
    im = ax.imshow(df.to_numpy().T / df.to_numpy().mean(axis=1)[:, None].T, cmap="Reds")
    fig.colorbar(im)

    ax.set_yticks(np.arange(len(design_names)), labels=design_names)
    ax.set_xticks(np.arange(len(criteria_names)), labels=criteria_names)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=20, ha="right")

    for i in range(len(design_names)):
        for j, criterion_name in enumerate(criteria_names):
            text = ax.text(
                j,
                i,
                f"{dfs[i].loc[criterion_name]['expected']:.1f}",
                ha="center",
                va="center",
                weight="bold" if criterion_name in selected_criteria else "normal",
            )

    save_plot(".", "tradeoff_values", type="png")
    plt.show()

    df["std_pct"] = df.std(axis=1) / df.mean(axis=1) * 100

    with pd.option_context(
        "display.max_rows", None, "display.max_columns", None, "display.width", None
    ):
        print(df)
