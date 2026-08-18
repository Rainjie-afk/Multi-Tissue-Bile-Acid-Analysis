

import argparse
import json
import os
import random
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.pipeline import Pipeline
from sklearn.model_selection import (
    StratifiedKFold,
    GridSearchCV,
)
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
    confusion_matrix,
    classification_report,
    roc_curve,
    auc,
)

try:
    from xgboost import XGBClassifier
except ImportError as exc:
    raise ImportError(
        "xgboost is required. Install it with: pip install xgboost"
    ) from exc





def set_seed(seed: int) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)





class RelativeAbundanceFilter(BaseEstimator, TransformerMixin):
\
\
\
\
\
\
\
\
\
\
\
\


    def __init__(
        self,
        min_prevalence=0.20,
        min_mean_abundance_percent=0.01,
    ):
        self.min_prevalence = min_prevalence
        self.min_mean_abundance_percent = min_mean_abundance_percent

    def fit(self, X, y=None):
        X_df = self._to_dataframe(X)

        prevalence = (X_df > 0).mean(axis=0)
        mean_abundance = X_df.mean(axis=0)

        support = (
            (prevalence >= self.min_prevalence)
            & (mean_abundance >= self.min_mean_abundance_percent)
        )


        if support.sum() == 0:
            warnings.warn(
                "No features passed the abundance filter. "
                "Keeping the single feature with the highest mean abundance."
            )
            support.iloc[np.argmax(mean_abundance.to_numpy())] = True

        self.feature_names_in_ = np.asarray(X_df.columns, dtype=object)
        self.support_mask_ = support.to_numpy(dtype=bool)
        self.selected_features_ = self.feature_names_in_[self.support_mask_]
        self.n_features_in_ = X_df.shape[1]
        self.n_features_out_ = int(self.support_mask_.sum())
        self.prevalence_ = prevalence.to_numpy()
        self.mean_abundance_percent_ = mean_abundance.to_numpy()

        return self

    def transform(self, X):
        X_df = self._to_dataframe(X)
        arr = X_df.to_numpy(dtype=float)
        return arr[:, self.support_mask_]

    def get_feature_names_out(self, input_features=None):
        return np.asarray(self.selected_features_, dtype=object)

    @staticmethod
    def _to_dataframe(X):
        if isinstance(X, pd.DataFrame):
            return X
        return pd.DataFrame(X)


class ArcsineSqrtPercentTransformer(BaseEstimator, TransformerMixin):
\
\
\
\
\


    def fit(self, X, y=None):
        return self

    def transform(self, X):
        X = np.asarray(X, dtype=float)

        if not np.isfinite(X).all():
            raise ValueError("Non-finite abundance values found.")

        if np.any(X < -1e-12) or np.any(X > 100 + 1e-12):
            raise ValueError(
                "ArcsineSqrtPercentTransformer expects percentages in [0, 100]."
            )

        proportion = np.clip(X / 100.0, 0.0, 1.0)
        return np.arcsin(np.sqrt(proportion))





def load_metadata(path, sample_col="Sample", group_col="Group"):
    meta = pd.read_csv(path, dtype={sample_col: str, group_col: str})

    required = {sample_col, group_col}
    missing = required.difference(meta.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")

    meta = meta[[sample_col, group_col]].copy()
    meta[sample_col] = meta[sample_col].astype(str).str.strip()
    meta[group_col] = meta[group_col].astype(str).str.strip()

    if meta[sample_col].duplicated().any():
        dup = meta.loc[meta[sample_col].duplicated(), sample_col].tolist()
        raise ValueError(f"Duplicate sample IDs in metadata: {dup}")

    return meta


def load_abundance_table(
    path,
    metadata,
    sample_col="Sample",
    feature_col="OUT",
    input_scale="counts",
):
\
\
\
\
\
\
\

    table = pd.read_csv(path, dtype={feature_col: str})

    if feature_col not in table.columns:
        raise ValueError(
            f"Feature ID column '{feature_col}' not found. "
            f"Available columns: {table.columns.tolist()[:10]} ..."
        )

    if table[feature_col].duplicated().any():
        raise ValueError("Feature IDs are not unique.")

    wanted_samples = metadata[sample_col].astype(str).tolist()
    abundance_columns = {str(c): c for c in table.columns if c != feature_col}

    missing_samples = [s for s in wanted_samples if s not in abundance_columns]
    if missing_samples:
        raise ValueError(
            "These metadata samples are absent from the abundance table: "
            + ", ".join(missing_samples)
        )


    selected_original_cols = [abundance_columns[s] for s in wanted_samples]

    numeric = table[selected_original_cols].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        bad = int(numeric.isna().sum().sum())
        raise ValueError(f"Abundance table contains {bad} non-numeric/NA values.")

    if (numeric < 0).any().any():
        raise ValueError("Negative abundances are not allowed.")

    feature_names = table[feature_col].astype(str).tolist()


    X = numeric.T
    X.index = wanted_samples
    X.columns = feature_names
    X.index.name = sample_col

    input_scale = input_scale.lower()

    if input_scale == "counts":
        totals = X.sum(axis=1)
        if (totals <= 0).any():
            bad_samples = totals[totals <= 0].index.tolist()
            raise ValueError(
                f"Samples with non-positive total counts: {bad_samples}"
            )
        X_percent = X.div(totals, axis=0) * 100.0

    elif input_scale == "percent":
        X_percent = X.astype(float)
        if (X_percent > 100 + 1e-8).any().any():
            raise ValueError(
                "Values >100 found although --input-scale percent was selected."
            )

    elif input_scale == "proportion":
        if (X > 1 + 1e-8).any().any():
            raise ValueError(
                "Values >1 found although --input-scale proportion was selected."
            )
        X_percent = X.astype(float) * 100.0

    else:
        raise ValueError(
            "input_scale must be one of: counts, percent, proportion"
        )

    return X_percent





def make_models(random_state, n_classes):
\
\
\


    common_preprocessing = [
        (
            "abundance_filter",
            RelativeAbundanceFilter(
                min_prevalence=0.20,
                min_mean_abundance_percent=0.001,
            ),
        ),
        ("arcsine_sqrt", ArcsineSqrtPercentTransformer()),
    ]

    models = {}


    lr = Pipeline(
        common_preprocessing
        + [
            ("scale", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    penalty="elasticnet",
                    solver="saga",
                    max_iter=5000,
                    class_weight="balanced",
                    random_state=random_state,
                ),
            ),
        ]
    )
    lr_grid = {
        "clf__C": [0.01, 0.1, 1.0],
        "clf__l1_ratio": [0.0, 0.5, 1.0],
    }
    models["ElasticNet_Logistic"] = (lr, lr_grid)


    rf = Pipeline(
        common_preprocessing
        + [
            (
                "clf",
                RandomForestClassifier(
                    n_estimators=300,
                    class_weight="balanced",
                    random_state=random_state,
                    n_jobs=1,
                ),
            )
        ]
    )
    rf_grid = {
        "clf__max_depth": [None, 3],
        "clf__min_samples_leaf": [1, 2],
        "clf__max_features": ["sqrt", 0.30],
    }
    models["Random_Forest"] = (rf, rf_grid)


    xgb = Pipeline(
        common_preprocessing
        + [
            (
                "clf",
                XGBClassifier(
                    objective="multi:softprob",
                    num_class=n_classes,
                    eval_metric="mlogloss",
                    n_estimators=150,
                    tree_method="hist",
                    random_state=random_state,
                    n_jobs=1,
                    verbosity=0,
                ),
            )
        ]
    )
    xgb_grid = {
        "clf__max_depth": [1, 2],
        "clf__learning_rate": [0.03, 0.10],
        "clf__subsample": [0.8],
        "clf__colsample_bytree": [0.5, 1.0],
    }
    models["XGBoost"] = (xgb, xgb_grid)




    dnn = Pipeline(
        common_preprocessing
        + [
            ("scale", StandardScaler()),
            (
                "clf",
                MLPClassifier(
                    hidden_layer_sizes=(16, 8),
                    activation="relu",
                    solver="lbfgs",
                    max_iter=2000,
                    random_state=random_state,
                ),
            ),
        ]
    )
    dnn_grid = {
        "clf__hidden_layer_sizes": [(8, 4), (16, 8)],
        "clf__alpha": [0.01, 0.1],
    }
    models["DNN_exploratory"] = (dnn, dnn_grid)

    return models





def extract_selected_features(best_pipeline):
    filt = best_pipeline.named_steps["abundance_filter"]
    return list(filt.get_feature_names_out())


def extract_model_importance(best_pipeline, model_name):
\
\
\
\
\
\

    feature_names = extract_selected_features(best_pipeline)
    clf = best_pipeline.named_steps["clf"]

    if model_name == "ElasticNet_Logistic":
        importance = np.mean(np.abs(clf.coef_), axis=0)

    elif model_name in {"Random_Forest", "XGBoost"}:
        importance = np.asarray(clf.feature_importances_, dtype=float)

    else:
        return None

    if len(importance) != len(feature_names):
        return None

    out = pd.DataFrame(
        {
            "Feature": feature_names,
            "Importance": importance,
        }
    )
    return out.sort_values("Importance", ascending=False, ignore_index=True)





def evaluate_models(
    X_percent,
    y,
    sample_ids,
    class_names,
    output_dir,
    n_splits=5,
    n_repeats=20,
    inner_splits=3,
    seed=0,
    n_jobs=-1,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_classes = len(class_names)
    all_metrics = []
    all_oof = []
    all_confusions = []
    all_selected = []
    all_importances = []
    all_best_params = []

    for repeat in range(n_repeats):
        repeat_seed = seed + repeat
        outer_cv = StratifiedKFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=repeat_seed,
        )

        print("\n" + "=" * 78)
        print(f"OUTER REPEAT {repeat + 1}/{n_repeats}  seed={repeat_seed}")
        print("=" * 78)


        model_defs = make_models(repeat_seed, n_classes)

        for model_name, (pipeline, param_grid) in model_defs.items():
            print(f"\nModel: {model_name}")

            oof_prob = np.full((len(y), n_classes), np.nan, dtype=float)
            oof_pred = np.full(len(y), -1, dtype=int)

            for fold, (train_idx, test_idx) in enumerate(
                outer_cv.split(X_percent, y), start=1
            ):
                X_train = X_percent.iloc[train_idx].copy()
                X_test = X_percent.iloc[test_idx].copy()
                y_train = y[train_idx]
                y_test = y[test_idx]


                inner_cv = StratifiedKFold(
                    n_splits=inner_splits,
                    shuffle=True,
                    random_state=repeat_seed * 100 + fold,
                )

                search = GridSearchCV(
                    estimator=pipeline,
                    param_grid=param_grid,
                    scoring="balanced_accuracy",
                    cv=inner_cv,
                    n_jobs=n_jobs,
                    refit=True,
                    error_score="raise",
                    return_train_score=False,
                )

                search.fit(X_train, y_train)
                best_model = search.best_estimator_

                prob = best_model.predict_proba(X_test)
                pred = np.argmax(prob, axis=1)

                oof_prob[test_idx, :] = prob
                oof_pred[test_idx] = pred

                selected_features = extract_selected_features(best_model)

                for feature in selected_features:
                    all_selected.append(
                        {
                            "Model": model_name,
                            "Repeat": repeat + 1,
                            "Fold": fold,
                            "Feature": feature,
                        }
                    )

                imp = extract_model_importance(best_model, model_name)
                if imp is not None:
                    imp = imp.copy()
                    imp["Model"] = model_name
                    imp["Repeat"] = repeat + 1
                    imp["Fold"] = fold
                    all_importances.append(imp)

                params_record = {
                    "Model": model_name,
                    "Repeat": repeat + 1,
                    "Fold": fold,
                    "Inner_best_score": search.best_score_,
                    "N_train": len(train_idx),
                    "N_test": len(test_idx),
                    "N_features_after_filter": len(selected_features),
                    "Best_params": json.dumps(search.best_params_, sort_keys=True),
                }
                all_best_params.append(params_record)

                print(
                    f"  Fold {fold}: "
                    f"train={len(train_idx)}, test={len(test_idx)}, "
                    f"features={len(selected_features)}, "
                    f"inner balanced accuracy={search.best_score_:.3f}"
                )

            if np.isnan(oof_prob).any() or np.any(oof_pred < 0):
                raise RuntimeError("Incomplete OOF predictions detected.")


            acc = accuracy_score(y, oof_pred)
            bacc = balanced_accuracy_score(y, oof_pred)
            macro_f1 = f1_score(y, oof_pred, average="macro")

            try:
                macro_auc = roc_auc_score(
                    y,
                    oof_prob,
                    multi_class="ovr",
                    average="macro",
                )
            except ValueError:
                macro_auc = np.nan

            all_metrics.append(
                {
                    "Model": model_name,
                    "Repeat": repeat + 1,
                    "Accuracy": acc,
                    "Balanced_accuracy": bacc,
                    "Macro_F1": macro_f1,
                    "Macro_OVR_AUC": macro_auc,
                }
            )

            cm = confusion_matrix(y, oof_pred, labels=np.arange(n_classes))
            for i, true_name in enumerate(class_names):
                for j, pred_name in enumerate(class_names):
                    all_confusions.append(
                        {
                            "Model": model_name,
                            "Repeat": repeat + 1,
                            "True_class": true_name,
                            "Predicted_class": pred_name,
                            "N": int(cm[i, j]),
                        }
                    )

            for i, sid in enumerate(sample_ids):
                row = {
                    "Model": model_name,
                    "Repeat": repeat + 1,
                    "Sample": sid,
                    "True_class": class_names[y[i]],
                    "Predicted_class": class_names[oof_pred[i]],
                }
                for k, cname in enumerate(class_names):
                    row[f"Prob_{cname}"] = oof_prob[i, k]
                all_oof.append(row)

            print(
                f"  OOF metrics: accuracy={acc:.3f}, "
                f"balanced_accuracy={bacc:.3f}, "
                f"macro_F1={macro_f1:.3f}, "
                f"macro_OVR_AUC={macro_auc:.3f}"
            )

    metrics_df = pd.DataFrame(all_metrics)
    oof_df = pd.DataFrame(all_oof)
    confusion_df = pd.DataFrame(all_confusions)
    selected_df = pd.DataFrame(all_selected)
    params_df = pd.DataFrame(all_best_params)

    metrics_df.to_csv(output_dir / "repeat_metrics.csv", index=False)
    oof_df.to_csv(output_dir / "oof_predictions.csv", index=False)
    confusion_df.to_csv(output_dir / "confusion_counts.csv", index=False)
    selected_df.to_csv(output_dir / "abundance_filter_selected_features.csv", index=False)
    params_df.to_csv(output_dir / "best_hyperparameters_by_fold.csv", index=False)

    if all_importances:
        importance_df = pd.concat(all_importances, ignore_index=True)
        importance_df.to_csv(
            output_dir / "outer_training_feature_importance.csv",
            index=False,
        )
    else:
        importance_df = pd.DataFrame()

    summary = (
        metrics_df.groupby("Model")
        .agg(
            Accuracy_mean=("Accuracy", "mean"),
            Accuracy_sd=("Accuracy", "std"),
            Balanced_accuracy_mean=("Balanced_accuracy", "mean"),
            Balanced_accuracy_sd=("Balanced_accuracy", "std"),
            Macro_F1_mean=("Macro_F1", "mean"),
            Macro_F1_sd=("Macro_F1", "std"),
            Macro_OVR_AUC_mean=("Macro_OVR_AUC", "mean"),
            Macro_OVR_AUC_sd=("Macro_OVR_AUC", "std"),
        )
        .reset_index()
    )
    summary.to_csv(output_dir / "model_performance_summary.csv", index=False)



    n_outer_fits = n_repeats * n_splits
    selection_frequency = (
        selected_df.groupby(["Model", "Feature"])
        .size()
        .rename("N_selected")
        .reset_index()
    )
    selection_frequency["Selection_frequency"] = (
        selection_frequency["N_selected"] / n_outer_fits
    )
    selection_frequency = selection_frequency.sort_values(
        ["Model", "Selection_frequency", "Feature"],
        ascending=[True, False, True],
    )
    selection_frequency.to_csv(
        output_dir / "abundance_filter_selection_frequency.csv",
        index=False,
    )


    if not importance_df.empty:
        importance_summary = (
            importance_df.groupby(["Model", "Feature"])
            .agg(
                Mean_importance=("Importance", "mean"),
                SD_importance=("Importance", "std"),
                N_outer_folds=("Importance", "size"),
            )
            .reset_index()
            .sort_values(["Model", "Mean_importance"], ascending=[True, False])
        )
        importance_summary.to_csv(
            output_dir / "feature_importance_summary.csv",
            index=False,
        )

    save_performance_plot(metrics_df, output_dir)
    save_mean_confusion_plots(confusion_df, class_names, output_dir)
    save_mean_roc_plots(oof_df, class_names, output_dir)

    return summary





def save_performance_plot(metrics_df, output_dir):
    metric_cols = [
        "Accuracy",
        "Balanced_accuracy",
        "Macro_F1",
        "Macro_OVR_AUC",
    ]

    means = metrics_df.groupby("Model")[metric_cols].mean()
    sds = metrics_df.groupby("Model")[metric_cols].std()

    x = np.arange(len(means.index))
    width = 0.18

    plt.figure(figsize=(10, 6))
    for i, metric in enumerate(metric_cols):
        plt.bar(
            x + (i - 1.5) * width,
            means[metric].to_numpy(),
            width=width,
            yerr=sds[metric].to_numpy(),
            capsize=3,
            label=metric,
        )

    plt.axhline(1 / 3, linestyle="--", linewidth=1, label="Chance accuracy (1/3)")
    plt.xticks(x, means.index, rotation=20, ha="right")
    plt.ylim(0, 1.05)
    plt.ylabel("Repeated-CV performance")
    plt.title("Three-class model comparison")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(output_dir / "model_performance_summary.pdf")
    plt.savefig(output_dir / "model_performance_summary.png", dpi=300)
    plt.close()


def save_mean_confusion_plots(confusion_df, class_names, output_dir):
    for model_name in confusion_df["Model"].unique():
        sub = confusion_df[confusion_df["Model"] == model_name]
        mean_cm = (
            sub.groupby(["True_class", "Predicted_class"])["N"]
            .mean()
            .unstack(fill_value=0)
            .reindex(index=class_names, columns=class_names)
        )


        row_sums = mean_cm.sum(axis=1).replace(0, np.nan)
        normalized = mean_cm.div(row_sums, axis=0)

        plt.figure(figsize=(5.5, 5))
        plt.imshow(normalized.to_numpy(), vmin=0, vmax=1, aspect="equal")
        plt.colorbar(label="Mean row-normalized proportion")
        plt.xticks(np.arange(len(class_names)), class_names, rotation=30)
        plt.yticks(np.arange(len(class_names)), class_names)
        plt.xlabel("Predicted class")
        plt.ylabel("True class")
        plt.title(f"Mean confusion matrix: {model_name}")

        for i in range(len(class_names)):
            for j in range(len(class_names)):
                val = normalized.iloc[i, j]
                plt.text(j, i, f"{val:.2f}", ha="center", va="center")

        plt.tight_layout()
        safe_name = model_name.replace(" ", "_")

        plt.savefig(
            output_dir / f"confusion_{safe_name}.pdf",
            bbox_inches="tight"
        )


        plt.savefig(
            output_dir / f"confusion_{safe_name}.tiff",
            dpi=800,
            bbox_inches="tight"
        )

plt.close()


def save_mean_roc_plots(oof_df, class_names, output_dir):
\
\
\

    prob_cols = [f"Prob_{c}" for c in class_names]

    for model_name in oof_df["Model"].unique():
        sub = oof_df[oof_df["Model"] == model_name].copy()

        mean_prob = (
            sub.groupby(["Sample", "True_class"], as_index=False)[prob_cols]
            .mean()
        )

        true_index = mean_prob["True_class"].map(
            {c: i for i, c in enumerate(class_names)}
        ).to_numpy()

        y_bin = label_binarize(
            true_index,
            classes=np.arange(len(class_names)),
        )

        plt.figure(figsize=(6.5, 5.5))

        aucs = []
        for k, cname in enumerate(class_names):
            fpr, tpr, _ = roc_curve(
                y_bin[:, k],
                mean_prob[f"Prob_{cname}"].to_numpy(),
            )
            class_auc = auc(fpr, tpr)
            aucs.append(class_auc)
            plt.plot(fpr, tpr, linewidth=2, label=f"{cname} (AUC={class_auc:.2f})")

        plt.plot([0, 1], [0, 1], "k--", linewidth=1)
        plt.xlim(0, 1)
        plt.ylim(0, 1.02)
        plt.xlabel("False positive rate")
        plt.ylabel("True positive rate")
        plt.title(
            f"One-vs-rest ROC: {model_name}\n"
            f"mean class AUC={np.mean(aucs):.2f}"
        )
        plt.legend(loc="lower right", fontsize=8)
        plt.tight_layout()

        safe_name = model_name.replace(" ", "_")
        plt.savefig(output_dir / f"ROC_{safe_name}.pdf")
        plt.close()





def parse_args():
    parser = argparse.ArgumentParser(
        description="Leakage-safe 3-class ML for NAFLD microbiome abundances."
    )

    parser.add_argument(
        "--metadata",
        default="NAFLD_Group_MCD2.csv",
        help="Metadata CSV containing Sample and Group columns.",
    )
    parser.add_argument(
        "--abundance",
        default="NAFLD_OTU(1).csv",
        help="Feature x sample abundance CSV.",
    )
    parser.add_argument(
        "--output",
        default="NAFLD_ML_results",
        help="Output directory.",
    )
    parser.add_argument("--sample-col", default="Sample")
    parser.add_argument("--group-col", default="Group")
    parser.add_argument("--feature-col", default="OUT")

    parser.add_argument(
        "--input-scale",
        choices=["counts", "percent", "proportion"],
        default="counts",
        help=(
            "Scale of the abundance input. The supplied NAFLD_OTU(1).csv "
            "contains counts, so 'counts' is the correct default."
        ),
    )

    parser.add_argument(
        "--min-prevalence",
        type=float,
        default=0.20,
        help="Minimum nonzero prevalence for training-fold feature filtering.",
    )
    parser.add_argument(
        "--min-mean-abundance-percent",
        type=float,
        default=0.0001,
        help="Minimum mean relative abundance percentage in a training fold.",
    )

    parser.add_argument(
        "--n-splits",
        type=int,
        default=5,
        help="Number of stratified outer folds. Default=5.",
    )
    parser.add_argument(
        "--n-repeats",
        type=int,
        default=20,
        help="Number of repeated outer CV runs. Default=20.",
    )
    parser.add_argument(
        "--inner-splits",
        type=int,
        default=3,
        help="Inner stratified folds for hyperparameter tuning. Default=3.",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="Parallel jobs used by GridSearchCV.",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    set_seed(args.seed)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    meta = load_metadata(
        args.metadata,
        sample_col=args.sample_col,
        group_col=args.group_col,
    )

    group_counts = meta[args.group_col].value_counts()
    if len(group_counts) != 3:
        warnings.warn(
            f"Expected 3 groups but found {len(group_counts)}: "
            f"{group_counts.to_dict()}"
        )

    if group_counts.min() < args.n_splits:
        raise ValueError(
            f"n_splits={args.n_splits} is too large. "
            f"The smallest class has only {group_counts.min()} samples."
        )

    if group_counts.min() - int(np.ceil(group_counts.min() / args.n_splits)) < args.inner_splits:
        warnings.warn(
            "The requested inner CV may be aggressive for the smallest "
            "outer-training class. Consider fewer inner folds if errors occur."
        )

    X_percent = load_abundance_table(
        args.abundance,
        metadata=meta,
        sample_col=args.sample_col,
        feature_col=args.feature_col,
        input_scale=args.input_scale,
    )






    global USER_MIN_PREVALENCE
    global USER_MIN_MEAN_ABUNDANCE_PERCENT
    USER_MIN_PREVALENCE = args.min_prevalence
    USER_MIN_MEAN_ABUNDANCE_PERCENT = args.min_mean_abundance_percent


    original_make_models = make_models

    def configured_make_models(random_state, n_classes):
        model_dict = original_make_models(random_state, n_classes)
        for _, (pipe, _) in model_dict.items():
            pipe.set_params(
                abundance_filter__min_prevalence=USER_MIN_PREVALENCE,
                abundance_filter__min_mean_abundance_percent=USER_MIN_MEAN_ABUNDANCE_PERCENT,
            )
        return model_dict

    globals()["make_models"] = configured_make_models


    class_names = sorted(meta[args.group_col].unique().tolist())
    class_to_int = {c: i for i, c in enumerate(class_names)}
    y = meta[args.group_col].map(class_to_int).to_numpy(dtype=int)
    sample_ids = meta[args.sample_col].astype(str).to_numpy()


    if list(X_percent.index.astype(str)) != list(sample_ids):
        raise RuntimeError("Sample order mismatch after loading.")

    if X_percent.shape[0] != len(meta):
        raise RuntimeError("Number of abundance samples does not match metadata.")


    X_percent.to_csv(output_dir / "relative_abundance_percent_all_samples.csv")




    full_arcsine_qc = np.arcsin(
        np.sqrt(np.clip(X_percent.to_numpy(dtype=float) / 100.0, 0.0, 1.0))
    )
    pd.DataFrame(
        full_arcsine_qc,
        index=X_percent.index,
        columns=X_percent.columns,
    ).to_csv(output_dir / "arcsine_sqrt_percent_div100_QC_only.csv")

    with open(output_dir / "class_mapping.json", "w", encoding="utf-8") as f:
        json.dump(class_to_int, f, indent=2)

    run_info = {
        "n_samples": int(X_percent.shape[0]),
        "n_features_before_filtering": int(X_percent.shape[1]),
        "group_counts": {str(k): int(v) for k, v in group_counts.items()},
        "input_scale": args.input_scale,
        "relative_abundance_unit_used_for_filtering": "percent",
        "arcsine_sqrt_formula": "arcsin(sqrt(relative_abundance_percent / 100))",
        "min_prevalence": args.min_prevalence,
        "min_mean_abundance_percent": args.min_mean_abundance_percent,
        "outer_cv": f"Repeated StratifiedKFold: {args.n_splits} folds x {args.n_repeats} repeats",
        "inner_cv": f"StratifiedKFold: {args.inner_splits} folds",
        "models": [
            "ElasticNet_Logistic",
            "Random_Forest",
            "XGBoost",
            "DNN_exploratory (2-hidden-layer MLP)",
        ],
        "note": (
            "All abundance filtering and arcsine-square-root transformation "
            "used for model fitting occur inside the CV pipeline."
        ),
    }

    with open(output_dir / "run_info.json", "w", encoding="utf-8") as f:
        json.dump(run_info, f, indent=2)

    print("\nDATA SUMMARY")
    print("-" * 60)
    print(f"Samples: {X_percent.shape[0]}")
    print(f"Features: {X_percent.shape[1]}")
    print("Groups:")
    print(group_counts.to_string())
    print(
        "\nPreprocessing used in every CV training pipeline:\n"
        f"  prevalence >= {args.min_prevalence}\n"
        f"  mean relative abundance >= {args.min_mean_abundance_percent}%\n"
        "  arcsin(sqrt(relative_abundance_percent / 100))"
    )
    print(
        f"\nOuter CV: stratified {args.n_splits}-fold x "
        f"{args.n_repeats} repeats"
    )
    print(f"Inner CV: stratified {args.inner_splits}-fold")

    summary = evaluate_models(
        X_percent=X_percent,
        y=y,
        sample_ids=sample_ids,
        class_names=class_names,
        output_dir=output_dir,
        n_splits=args.n_splits,
        n_repeats=args.n_repeats,
        inner_splits=args.inner_splits,
        seed=args.seed,
        n_jobs=args.n_jobs,
    )

    print("\nFINAL REPEATED-CV SUMMARY")
    print("=" * 78)
    print(summary.to_string(index=False))
    print(f"\nResults written to: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
