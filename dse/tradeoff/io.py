import numpy as np
import pandas as pd


def load_sheets(file, design_names, selected_only=True):
    # Load scores
    df_scoring = pd.read_excel(file, sheet_name="Scoring")
    df_scoring = df_scoring.rename(
        columns={"Category": "category", "Numerical score": "score"}
    ).set_index("category", drop=True)[["score"]]

    # Load weights
    df_weights = pd.read_excel(file, sheet_name="Criteria")

    # Remove rows that are not criteria
    first_empty_row_idx = df_weights["Criterion"].isnull().idxmax()
    df_weights = df_weights.iloc[:first_empty_row_idx]

    df_weights = df_weights.rename(
        columns={"Criterion": "criterion", "Weight": "weight", "Selected": "selected"}
    ).set_index("criterion", drop=True)
    df_weights["weight"] = df_weights["weight"].round(1)
    selected_criteria = df_weights[df_weights["selected"] == "x"].index

    # Extract score ranges per criterion
    score_categories = {
        criterion: {category: [None, None] for category in df_scoring.index}
        for criterion in selected_criteria
    }
    for criterion, categories in score_categories.items():
        for category in list(categories.keys())[:-1]:
            row = df_weights.loc[criterion]
            boundary = row.iloc[np.where(df_weights.columns == category)[0][0] - 1]
            score_categories[criterion][category][0] = boundary

        for left, right in zip(list(categories.keys())[1:], list(categories.keys())[:-1]):
            score_categories[criterion][left][1] = score_categories[criterion][right][0]

    df_weights = df_weights[["weight", "selected"]]

    # Load criteria values for each design
    dfs = []
    for design_name in design_names:
        df = pd.read_excel(file, sheet_name=design_name)
        df = df.drop("Notes", axis=1)
        df = df.iloc[: len(df_weights.index)]
        df = df.rename(
            columns={
                "Criterion": "criterion",
                "Expected": "expected",
                "Worst": "worst",
                "Best": "best",
            }
        ).set_index("criterion", drop=True)

        if selected_only:
            df = df.loc[selected_criteria]

        dfs.append(df)

    if selected_only:
        df_weights = df_weights.loc[selected_criteria]

    return dfs, df_weights, df_scoring, score_categories
