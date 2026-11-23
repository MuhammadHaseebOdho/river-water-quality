# %% Imports
import os
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import scikit_posthocs as sp

os.makedirs("results", exist_ok=True)


# %% Load and inspect the data
df = pd.read_csv("waterquality.csv", sep=";")
variables = ["nitrate_mgL", "oxygen_mgL"]
seasons = ["spring", "summer", "autumn", "winter"]
sections = sorted(df["river_section"].unique())
df["saison"] = pd.Categorical(df["saison"], categories=seasons, ordered=True)
df["river_section"] = pd.Categorical(df["river_section"], categories=sections, ordered=True)

print(df.head())
print("\nData types:\n", df.dtypes)
print("\nShape:", df.shape)
print("\nMissing values:\n", df.isna().sum())
print("Duplicate rows:", df.duplicated().sum())
print("Duplicate IDs:", df["ID"].duplicated().sum())
counts = pd.crosstab(df["river_section"], df["saison"])
print("\nMeasurements by section and season:")
print(counts)
counts.to_csv("results/sample_counts.csv")

# %% Descriptive analysis
summary_names = {"count": "n", "std": "sd", "25%": "Q1", "50%": "median", "75%": "Q3"}
overall = df[variables].describe().T.rename(columns=summary_names)
print("\nOverall summary:")
print(overall)
overall.to_csv("results/descriptive_overall.csv")

by_season = df.groupby("saison", observed=True)[variables].describe()
by_season = by_season.rename(columns=summary_names, level=1)
print("\nSummary by season:")
print(by_season)
by_season.to_csv("results/descriptive_by_season.csv")

by_section = df.groupby("river_section", observed=True)[variables].describe()
by_section = by_section.rename(columns=summary_names, level=1)
print("\nSummary by section:")
print(by_section)
by_section.to_csv("results/descriptive_by_section.csv")

means = df.groupby(["river_section", "saison"], observed=True)[variables].mean().unstack("saison")
print("\nMean concentrations by section and season:")
print(means)
means.to_csv("results/section_season_means.csv")

# %% Boxplots
sns.set_theme(style="whitegrid", context="paper", font_scale=1.6)
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
sns.boxplot(data=df, x="saison", y="nitrate_mgL", color="#83b8c9", ax=axes[0, 0])
sns.boxplot(data=df, x="river_section", y="nitrate_mgL", color="#83b8c9", ax=axes[0, 1])
sns.boxplot(data=df, x="saison", y="oxygen_mgL", color="#83b8c9", ax=axes[1, 0])
sns.boxplot(data=df, x="river_section", y="oxygen_mgL", color="#83b8c9", ax=axes[1, 1])

axes[0, 0].set(title="Nitrate by season", xlabel="Season", ylabel="Nitrate (mg/L)")
axes[0, 1].set(title="Nitrate by river section", xlabel="River section", ylabel="Nitrate (mg/L)")
axes[1, 0].set(title="Oxygen by season", xlabel="Season", ylabel="Oxygen (mg/L)")
axes[1, 1].set(title="Oxygen by river section", xlabel="River section", ylabel="Oxygen (mg/L)")

# A1 = headwaters, A2 = upper course, A3 = middle course, A4 = lower course, A5 = mouth
section_labels = ["A1", "A2", "A3", "A4", "A5"]
season_labels = ["Spring", "Summer", "Autumn", "Winter"]
axes[0, 0].set_xticks(range(4), season_labels)
axes[1, 0].set_xticks(range(4), season_labels)
axes[0, 1].set_xticks(range(5), section_labels)
axes[1, 1].set_xticks(range(5), section_labels)
fig.tight_layout()
fig.savefig("fig1_boxplots.png", dpi=300, bbox_inches="tight")

# %% Assumption checks
shapiro_rows = []
levene_rows = []
# RQ1: check spring and autumn separately for each concentration.
for variable in variables:
    spring = df.loc[df["saison"] == "spring", variable]
    autumn = df.loc[df["saison"] == "autumn", variable]
    spring_check = stats.shapiro(spring)
    autumn_check = stats.shapiro(autumn)
    shapiro_rows.append(["RQ1", variable, "spring", len(spring), spring_check.statistic, spring_check.pvalue])
    shapiro_rows.append(["RQ1", variable, "autumn", len(autumn), autumn_check.statistic, autumn_check.pvalue])
    variance_check = stats.levene(spring, autumn, center="median")
    levene_rows.append(["RQ1", variable, variance_check.statistic, variance_check.pvalue])

# RQ2: check all five sections for each concentration.
for variable in variables:
    groups = []
    for section in sections:
        values = df.loc[df["river_section"] == section, variable]
        groups.append(values)
        normality_check = stats.shapiro(values)
        shapiro_rows.append(["RQ2", variable, section, len(values), normality_check.statistic, normality_check.pvalue])
    # The asterisk passes the five groups as separate samples.
    variance_check = stats.levene(*groups, center="median")
    levene_rows.append(["RQ2", variable, variance_check.statistic, variance_check.pvalue])

shapiro = pd.DataFrame(shapiro_rows, columns=["question", "variable", "group", "n", "W", "p_value"])
levene = pd.DataFrame(levene_rows, columns=["question", "variable", "statistic", "p_value"])
print("\nShapiro-Wilk normality checks:")
print(shapiro)
shapiro.to_csv("results/shapiro_wilk.csv", index=False)
print("\nMedian-centred Levene variance checks:")
print(levene)
levene.to_csv("results/levene.csv", index=False)

# %% RQ1: Spring versus autumn
# All differences and effect sizes use spring minus autumn.
# The two primary tests use Bonferroni alpha = 0.05 / 2.
welch_rows = []
for variable in variables:
    spring = df.loc[df["saison"] == "spring", variable]
    autumn = df.loc[df["saison"] == "autumn", variable]
    result = stats.ttest_ind(spring, autumn, equal_var=False)
    ci = result.confidence_interval(confidence_level=0.95)
    n_spring = len(spring)
    n_autumn = len(autumn)
    # Cohen's d: mean difference divided by the pooled standard deviation.
    pooled_sd = np.sqrt(((n_spring - 1) * spring.var(ddof=1) + (n_autumn - 1) * autumn.var(ddof=1)) / (n_spring + n_autumn - 2))
    difference = spring.mean() - autumn.mean()
    cohens_d = difference / pooled_sd
    welch_rows.append({
        "variable": variable,
        "n_spring": n_spring,
        "n_autumn": n_autumn,
        "mean_spring": spring.mean(),
        "mean_autumn": autumn.mean(),
        "mean_difference": difference,
        "difference_ci95_low": ci.low,
        "difference_ci95_high": ci.high,
        "cohens_d": cohens_d,
        "t": result.statistic,
        "df": result.df,
        "p_value": result.pvalue,
        "p_bonferroni": min(2 * result.pvalue, 1.0),
        "reject_at_0.025": result.pvalue < 0.025,
    })

welch = pd.DataFrame(welch_rows)
print("\nSpring versus autumn: Welch t-tests")
print(welch)
welch.to_csv("results/rq1_welch.csv", index=False)

# %% RQ2: Nitrate by river section
nitrate_groups = []
for section in sections:
    values = df.loc[df["river_section"] == section, "nitrate_mgL"]
    nitrate_groups.append(values)
anova = stats.f_oneway(*nitrate_groups)
grand_mean = df["nitrate_mgL"].mean()
ss_between = 0
ss_within = 0
for group in nitrate_groups:
    ss_between += len(group) * (group.mean() - grand_mean) ** 2
    ss_within += ((group - group.mean()) ** 2).sum()
eta_squared = ss_between / (ss_between + ss_within)
k = len(sections)
n_total = len(df)
anova_table = pd.DataFrame({
    "source": ["river_section", "within", "total"],
    "sum_sq": [ss_between, ss_within, ss_between + ss_within],
    "df": [k - 1, n_total - k, n_total - 1],
    "mean_sq": [ss_between / (k - 1), ss_within / (n_total - k), np.nan],
    "F": [anova.statistic, np.nan, np.nan],
    "p_value": [anova.pvalue, np.nan, np.nan],
    "eta_squared": [eta_squared, np.nan, np.nan],
}).set_index("source")
print("\nNitrate: one-way ANOVA")
print(anova_table)
anova_table.to_csv("results/rq2_nitrate_anova.csv")

tukey = pairwise_tukeyhsd(df["nitrate_mgL"], df["river_section"], alpha=0.05)
pairs = list(combinations(tukey.groupsunique, 2))
tukey_table = pd.DataFrame({
    "group1": [pair[0] for pair in pairs], "group2": [pair[1] for pair in pairs],
    "mean_difference_group2_minus_group1": tukey.meandiffs,
    "ci95_low": tukey.confint[:, 0], "ci95_high": tukey.confint[:, 1],
    "p_tukey": tukey.pvalues, "reject_at_0.05": tukey.reject,
})
print("\nNitrate: Tukey comparisons")
print(tukey_table)
tukey_table.to_csv("results/rq2_nitrate_tukey.csv", index=False)

# %% RQ2: Oxygen by river section
oxygen_groups = []
for section in sections:
    values = df.loc[df["river_section"] == section, "oxygen_mgL"]
    oxygen_groups.append(values)
kruskal = stats.kruskal(*oxygen_groups)
kruskal_table = pd.DataFrame({"H": [kruskal.statistic], "df": [len(sections) - 1], "p_value": [kruskal.pvalue]})
print("\nOxygen: Kruskal-Wallis test")
print(kruskal_table)
kruskal_table.to_csv("results/rq2_oxygen_kruskal.csv", index=False)

dunn_raw = sp.posthoc_dunn(df, val_col="oxygen_mgL", group_col="river_section", p_adjust=None)
dunn_holm = sp.posthoc_dunn(df, val_col="oxygen_mgL", group_col="river_section", p_adjust="holm")
dunn_rows = []
for group1, group2 in combinations(sections, 2):
    p_value = dunn_raw.loc[group1, group2]
    p_holm = dunn_holm.loc[group1, group2]
    dunn_rows.append([group1, group2, p_value, p_holm, p_holm < 0.05])
dunn_table = pd.DataFrame(dunn_rows, columns=["group1", "group2", "p_value", "p_holm", "reject_at_0.05"])
print("\nOxygen: Dunn comparisons with Holm correction")
print(dunn_table)
dunn_table.to_csv("results/rq2_oxygen_dunn_holm.csv", index=False)

# %% RQ2: Two-way ANOVA with section and season
# Type II sums of squares; partial eta squared = SS_effect / (SS_effect + SS_residual).
effects = ["C(river_section)", "C(saison)", "C(river_section):C(saison)"]
effect_names = ["river_section", "saison", "river_section:saison"]
two_way_rows = []
two_way_checks = []
for variable in variables:
    model = smf.ols(f"{variable} ~ C(river_section) * C(saison)", data=df).fit()
    table = sm.stats.anova_lm(model, typ=2)
    if variable == "oxygen_mgL":
        oxygen_model = model
    ss_residual = table.loc["Residual", "sum_sq"]
    for effect, name in zip(effects, effect_names):
        ss_effect = table.loc[effect, "sum_sq"]
        partial_eta_squared = ss_effect / (ss_effect + ss_residual)
        two_way_rows.append([variable, name, table.loc[effect, "df"], model.df_resid,
                             table.loc[effect, "F"], table.loc[effect, "PR(>F)"], partial_eta_squared])
    # Residual normality and equal variances across the 20 section-season cells.
    cells = []
    for section in sections:
        for season in seasons:
            cells.append(df.loc[(df["river_section"] == section) & (df["saison"] == season), variable])
    normality_check = stats.shapiro(model.resid)
    variance_check = stats.levene(*cells, center="median")
    two_way_checks.append([variable, normality_check.statistic, normality_check.pvalue,
                           variance_check.statistic, variance_check.pvalue])

two_way = pd.DataFrame(two_way_rows, columns=["variable", "effect", "df", "df_error", "F", "p_value",
                                              "partial_eta_squared"])
print("\nTwo-way ANOVA (type II sums of squares):")
print(two_way)
two_way.to_csv("results/rq2_two_way_anova.csv", index=False)

two_way_check_table = pd.DataFrame(two_way_checks, columns=["variable", "shapiro_W", "shapiro_p",
                                                            "levene_statistic", "levene_p"])
print("\nTwo-way ANOVA residual checks:")
print(two_way_check_table)
two_way_check_table.to_csv("results/rq2_two_way_checks.csv", index=False)

# Oxygen: section effect within each season, using the two-way error term.
# Bonferroni alpha = 0.05 / 4 for the four seasons.
mse = oxygen_model.mse_resid
df_error = oxygen_model.df_resid
season_rows = []
for season in seasons:
    season_values = df.loc[df["saison"] == season, "oxygen_mgL"]
    cell_means = df[df["saison"] == season].groupby("river_section", observed=True)["oxygen_mgL"].mean()
    ss_section = 45 * ((cell_means - season_values.mean()) ** 2).sum()
    f_value = (ss_section / (k - 1)) / mse
    p_value = stats.f.sf(f_value, k - 1, df_error)
    season_rows.append([season, f_value, k - 1, df_error, p_value, p_value < 0.05 / 4])
section_by_season = pd.DataFrame(season_rows, columns=["season", "F", "df", "df_error", "p_value",
                                                       "reject_at_0.0125"])
print("\nOxygen: section effect within each season")
print(section_by_season)
section_by_season.to_csv("results/rq2_oxygen_section_by_season.csv", index=False)

# %% Nitrate ANOVA model checks
# In one-way ANOVA, each fitted value is its section's sample mean.
fitted = df.groupby("river_section", observed=True)["nitrate_mgL"].transform("mean")
residuals = df["nitrate_mgL"] - fitted

sns.set_theme(style="whitegrid", context="paper", font_scale=1.6)
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
axes[0].scatter(fitted, residuals, color="#397b91", alpha=0.4, s=15)
axes[0].axhline(0, color="black", linestyle="--", linewidth=1)
axes[0].set(title="(a) Residuals versus fitted values", xlabel="Fitted nitrate (mg/L)",
            ylabel="Residual (mg/L)")

# probplot returns the normal quantiles, the ordered residuals and a fitted reference line.
(quantiles, ordered), (slope, intercept, r) = stats.probplot(residuals, dist="norm")
axes[1].scatter(quantiles, ordered, color="#397b91", alpha=0.6, s=15)
axes[1].plot(quantiles, intercept + slope * quantiles, color="black", linewidth=1)
axes[1].set(title="(b) Normal Q–Q plot", xlabel="Theoretical normal quantiles",
            ylabel="Ordered residuals (mg/L)")
fig.tight_layout()
fig.savefig("fig2_nitrate_anova_diagnostics.png", dpi=300, bbox_inches="tight")
