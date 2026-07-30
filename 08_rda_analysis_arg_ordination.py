"""
Redundancy Analysis (RDA) Statistics and ARG Ordination Plot
Peshawar Wastewater Metagenomics -- Figure 8A (RDA biplot, ARG profiles)

Parses BioinCloud's RDA/DCA/PERMANOVA/envfit output tables, prints the summary
statistics reported in the Results (DCA gradient length, PERMANOVA variance
explained, envfit r^2/p for sector and temporal variables, sector centroids),
and reproduces the ARG RDA biplot (Figure 8A) with sector confidence ellipses,
top-loading ARG feature vectors, and envfit vectors.

The equivalent MGE RDA biplot (Figure 8B) uses the same BioinCloud RDA output
structure for the MGE ordination and was plotted the same way; only the ARG
panel's plotting code is included here.

Input (six tab-separated BioinCloud RDA/DCA output files):
    DCA_output.xls
    Samecity3targetspeshawar_RDA_envfit.xls
    Samecity3targetspeshawar_RDA_Factors_PERMANOVA.xls
    Samecity3targetspeshawar_RDA_features.xls
    Samecity3targetspeshawar_RDA_Group_PERMANOVA.xls
    Samecity3targetspeshawar_RDA_sample.xls

Output:
    Fig_RDA_ARG.svg / .png (300 dpi)
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse

DATA_DIR = "data"
OUT_DIR = "output"

# ARO -> gene name, for labelling the top-loading feature vectors
ARO_NAMES = {
    "ARO:3000491": "mef(A)", "ARO:3002983": "tet(X4)", "ARO:3001777": "VIM-1",
    "ARO:3000498": "dfr(A1)", "ARO:3000794": "aac(6')-Ib", "ARO:3004089": "mcr-3",
    "ARO:3002675": "APH(3')-Ia", "ARO:3003841": "OXA-427", "ARO:3000237": "TEM-15",
    "ARO:3000191": "SHV-1", "ARO:3000183": "CTX-M-3", "ARO:3000796": "aph(3'')-Ib",
    "ARO:3001214": "ANT(3'')-Ia", "ARO:3004054": "mcr-9", "ARO:3001329": "NDM-5",
}

SECTOR_COLORS = {"Hospital": "#D62728", "Community": "#1F77B4", "Slaughterhouse": "#FF7F0E"}
SECTOR_GROUPS = {
    "Hospital": ["PHW1", "PHW2"],
    "Community": ["PCW1", "PCW2"],
    "Slaughterhouse": ["PSLW1", "PSLW2"],
}
SAMPLE_COLORS = {s: SECTOR_COLORS[sec] for sec, ss in SECTOR_GROUPS.items() for s in ss}


def print_summary_stats():
    dca = pd.read_csv(f"{DATA_DIR}/DCA_output.xls", sep="\t", index_col=0)
    print("=== DCA GRADIENT LENGTHS ===")
    print("DCA1 range:", dca["DCA1"].max() - dca["DCA1"].min())
    print("DCA2 range:", dca["DCA2"].max() - dca["DCA2"].min())
    print("Interpretation: gradient length < 3 SD -> RDA appropriate; > 4 SD -> CCA preferred\n")

    perm_f = pd.read_csv(f"{DATA_DIR}/Samecity3targetspeshawar_RDA_Factors_PERMANOVA.xls", sep="\t")
    perm_g = pd.read_csv(f"{DATA_DIR}/Samecity3targetspeshawar_RDA_Group_PERMANOVA.xls", sep="\t")

    for label, perm in [("Factors (source + temporal)", perm_f), ("Group (sector only)", perm_g)]:
        total_var = perm["Variance"].sum()
        model_row = perm[perm["Unnamed: 0"] == "Model"]
        model_var = model_row["Variance"].values[0]
        print(f"=== PERMANOVA -- {label} ===")
        print(perm.to_string())
        print(f"Variance explained: {100 * model_var / total_var:.1f}%")
        print(f"F = {model_row['F'].values[0]:.3f}, p = {model_row['Pr(>F)'].values[0]:.3f}\n")

    envfit = pd.read_csv(f"{DATA_DIR}/Samecity3targetspeshawar_RDA_envfit.xls", sep="\t")
    print("=== ENVFIT RESULTS ===")
    print(envfit.to_string())
    source_row = envfit[envfit["Unnamed: 0"] == "source"]
    temporal_row = envfit[envfit["Unnamed: 0"] == "temporal"]
    print(f"\nSource: r2 = {source_row['r2'].values[0]:.3f}, p = {source_row['p value'].values[0]:.3f}")
    print(f"Temporal: r2 = {temporal_row['r2'].values[0]:.3f}, p = {temporal_row['p value'].values[0]:.3f}\n")

    samples = pd.read_csv(f"{DATA_DIR}/Samecity3targetspeshawar_RDA_sample.xls", sep="\t", index_col=0)
    print("=== SAMPLE RDA SCORES & SECTOR CENTROIDS ===")
    print(samples.to_string())
    for sector, members in SECTOR_GROUPS.items():
        c = samples.loc[members]
        print(f"{sector} centroid: RDA1={c['RDA1'].mean():.3f}, RDA2={c['RDA2'].mean():.3f}")

    features = pd.read_csv(f"{DATA_DIR}/Samecity3targetspeshawar_RDA_features.xls", sep="\t", index_col=0)
    features["magnitude"] = np.sqrt(features["RDA1"] ** 2 + features["RDA2"] ** 2)
    print(f"\n=== TOP 15 ARG FEATURE VECTORS (by magnitude, of {len(features)} total) ===")
    print(features.nlargest(15, "magnitude").to_string())


def plot_arg_rda_biplot():
    samples = pd.read_csv(f"{DATA_DIR}/Samecity3targetspeshawar_RDA_sample.xls", sep="\t", index_col=0)
    features = pd.read_csv(f"{DATA_DIR}/Samecity3targetspeshawar_RDA_features.xls", sep="\t", index_col=0)
    envfit = pd.read_csv(f"{DATA_DIR}/Samecity3targetspeshawar_RDA_envfit.xls", sep="\t", index_col=0)
    features["magnitude"] = np.sqrt(features["RDA1"] ** 2 + features["RDA2"] ** 2)

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 11, "axes.labelweight": "bold",
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 300, "savefig.dpi": 300,
    })

    fig, ax = plt.subplots(figsize=(7, 7))

    # Sector confidence ellipses (simple bounding circle around each sector's points)
    for sec, samps in SECTOR_GROUPS.items():
        pts = samples.loc[samps, ["RDA1", "RDA2"]].values
        cx, cy = pts.mean(axis=0)
        radius = max(np.sqrt((pts[:, 0] - cx) ** 2 + (pts[:, 1] - cy) ** 2)) + 1.2
        circle = plt.Circle((cx, cy), radius, color=SECTOR_COLORS[sec], alpha=0.08,
                             linewidth=1.2, linestyle="--", edgecolor=SECTOR_COLORS[sec],
                             fill=True, zorder=1)
        ax.add_patch(circle)

    # Top 12 ARG feature vectors
    top_feat = features.nlargest(12, "magnitude")
    scale = 18
    for aro, row in top_feat.iterrows():
        gene = ARO_NAMES.get(aro, aro.replace("ARO:", ""))
        x, y = row["RDA1"] * scale, row["RDA2"] * scale
        ax.annotate("", xy=(x, y), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color="#888888", lw=1.0,
                                    mutation_scale=8, alpha=0.7))
        ha = "left" if x >= 0 else "right"
        offset = 0.3 if x >= 0 else -0.3
        ax.text(x + offset, y, gene, fontsize=7, color="#444444",
                fontstyle="italic", va="center", ha=ha, zorder=6)

    # Envfit vectors (sector, temporal)
    env_scale = 6
    for var, row in envfit.iterrows():
        ex, ey = row["RDA1"] * env_scale, row["RDA2"] * env_scale
        ax.annotate("", xy=(ex, ey), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color="#2A2A2A", lw=1.8, mutation_scale=11))
        label = "Sector" if var == "source" else "Temporal"
        sig = "*" if row["p value"] < 0.05 else ""
        ax.text(ex * 1.12, ey * 1.12, f"{label}\n(r2={row['r2']:.2f}{sig})",
                fontsize=8, color="#2A2A2A", fontweight="bold", ha="center", va="center", zorder=7)

    # Sample points
    for samp, row in samples.iterrows():
        col = SAMPLE_COLORS[samp]
        ax.scatter(row["RDA1"], row["RDA2"], s=120, color=col, edgecolors="white",
                   linewidths=1.0, zorder=5)
        ax.text(row["RDA1"] + 0.2, row["RDA2"] + 0.2, samp, fontsize=8, color=col,
                fontstyle="italic", fontweight="bold", zorder=6)

    ax.axhline(0, color="#CCCCCC", lw=0.6, zorder=0)
    ax.axvline(0, color="#CCCCCC", lw=0.6, zorder=0)
    ax.set_xlabel("RDA1")
    ax.set_ylabel("RDA2")
    ax.set_aspect("equal")

    all_x = list(samples["RDA1"]) + [r * scale for r in top_feat["RDA1"]]
    all_y = list(samples["RDA2"]) + [r * scale for r in top_feat["RDA2"]]
    pad = 1.5
    ax.set_xlim(min(all_x) - pad, max(all_x) + pad)
    ax.set_ylim(min(all_y) - pad, max(all_y) + pad)

    ax.text(0.03, 0.98, "PERMANOVA: Sector F = 0.865, p = 0.733\nSector explains 36.6% of ARG variance",
            transform=ax.transAxes, fontsize=8, va="top", ha="left", color="#333333",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#CCCCCC", alpha=0.9))

    patches = [mpatches.Patch(color=c, label=s) for s, c in SECTOR_COLORS.items()]
    ax.legend(handles=patches, loc="lower right", framealpha=0.9, edgecolor="#CCCCCC",
              fontsize=9, title="Sector", title_fontsize=9)

    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    plt.savefig(f"{OUT_DIR}/Fig_RDA_ARG.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUT_DIR}/Fig_RDA_ARG.png", format="png", dpi=300, bbox_inches="tight")
    print(f"ARG RDA figure saved to {OUT_DIR}/")


if __name__ == "__main__":
    print_summary_stats()
    plot_arg_rda_biplot()
