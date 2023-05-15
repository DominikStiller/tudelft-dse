import pandas as pd


def load_sheets(file, design_names):
    # Load weights
    df_weights = pd.read_excel(file, sheet_name="Criteria")
    df_weights = df_weights.rename(
        columns={"Criterion": "criterion", "Weight": "weight", "Selected": "selected"}
    ).set_index("criterion", drop=True)[["weight", "selected"]]
    selected_criteria = df_weights[df_weights["selected"] == "x"].index
    df_weights = df_weights.loc[selected_criteria]

    # Load scores
    df_scoring = pd.read_excel(file, sheet_name="Scoring")
    df_scoring = df_scoring.rename(
        columns={"Category": "category", "Numerical score": "score"}
    ).set_index("category", drop=True)

    # Load criteria values for each design
    dfs = []
    for design_name in design_names:
        df = pd.read_excel(file, sheet_name=design_name)
        df = df.drop("Notes", axis=1)
        df = df.rename(
            columns={
                "Criterion": "criterion",
                "Expected": "expected",
                "Score for Expected": "expected_score",
                "Worst": "worst",
                "Score for Worst": "worst_score",
                "Best": "best",
                "Score for Best": "best_score",
            }
        ).set_index("criterion", drop=True)
        df = df.loc[selected_criteria]

        df["expected_score"] = df["expected_score"].apply(lambda s: df_scoring.loc[s]["score"])
        df["worst_score"] = df["worst_score"].apply(lambda s: df_scoring.loc[s]["score"])
        df["best_score"] = df["best_score"].apply(lambda s: df_scoring.loc[s]["score"])

        dfs.append(df)

    return dfs, df_weights
