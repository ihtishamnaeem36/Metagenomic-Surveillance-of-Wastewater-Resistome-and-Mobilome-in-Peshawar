"""
Procrustes Superimposition Plot
Peshawar Wastewater Metagenomics -- Figure 7

Plots the Procrustes superimposition of the ARG and MGE community ordinations
(target = ARG PCoA, rotated = MGE PCoA), with residual vectors connecting each
sample's position in the two ordinations, coloured by sector.

The Procrustes rotation itself (m^2, p-value, per-sample residuals) was
computed by the BioinCloud platform / vegan::protest() in R; this script only
reproduces the publication figure from that output.

Input:
    procrustes_result.txt   Tab-separated: samples, groups, X1, X2 (ARG
                            ordination coordinates), PC1, PC2 (rotated MGE
                            ordination coordinates)

Output:
    Fig_Procrustes.svg / .png (300 dpi)
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

DATA_FILE = "data/procrustes_result.txt"
OUT_DIR = "output"

COLORS = {
    "Hospital": "#D62728",
    "Community": "#1F77B4",
    "Slaughterhouse": "#FF7F0E",
}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica"],
    "font.size": 10, "axes.labelsize": 11, "axes.labelweight": "bold",
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300,
})


def main():
    df = pd.read_csv(DATA_FILE, sep="\t")

    arg = df[["X1", "X2"]].values      # ARG ordination (target) -- open circles
    mge = df[["PC1", "PC2"]].values    # MGE ordination (rotated) -- filled circles
    samples = df["samples"].tolist()
    groups = df["groups"].tolist()

    fig, ax = plt.subplots(figsize=(6.5, 6))

    for i, (s, g) in enumerate(zip(samples, groups)):
        col = COLORS[g]
        ax0, ay0 = arg[i]
        mx0, my0 = mge[i]

        ax.annotate("", xy=(mx0, my0), xytext=(ax0, ay0),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.4,
                                    mutation_scale=10, alpha=0.85))
        ax.scatter(ax0, ay0, s=90, facecolors="white", edgecolors=col,
                   linewidths=1.8, zorder=4)
        ax.scatter(mx0, my0, s=90, color=col, edgecolors="white",
                   linewidths=0.6, zorder=5)
        ax.text(mx0 + 0.025, my0 + 0.025, s, fontsize=7.5, color=col,
                fontstyle="italic", va="bottom", ha="left", zorder=6)

    ax.axhline(0, color="#CCCCCC", lw=0.6, zorder=0)
    ax.axvline(0, color="#CCCCCC", lw=0.6, zorder=0)
    ax.set_xlabel("Procrustes Dimension 1 (PCoA)")
    ax.set_ylabel("Procrustes Dimension 2 (PCoA)")
    ax.set_aspect("equal")

    all_pts = np.vstack([arg, mge])
    pad = 0.12
    ax.set_xlim(all_pts[:, 0].min() - pad, all_pts[:, 0].max() + pad)
    ax.set_ylim(all_pts[:, 1].min() - pad, all_pts[:, 1].max() + pad)

    ax.text(0.03, 0.97, r"$M^2$ = 0.744,  $p$ = 0.756 (999 permutations)",
            transform=ax.transAxes, fontsize=9, va="top", ha="left", color="#333333",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#CCCCCC", alpha=0.85))

    sector_patches = [mpatches.Patch(facecolor=COLORS[g], label=g)
                      for g in ["Hospital", "Community", "Slaughterhouse"]]
    sym_open = ax.scatter([], [], s=70, facecolors="white", edgecolors="#555555",
                          linewidths=1.5, label="ARG ordination (target)")
    sym_fill = ax.scatter([], [], s=70, color="#555555", edgecolors="white",
                          linewidths=0.6, label="MGE ordination (rotated)")
    ax.legend(handles=sector_patches + [sym_open, sym_fill], loc="lower right",
              framealpha=0.92, edgecolor="#CCCCCC", fontsize=8,
              title="Sector / Symbol", title_fontsize=8, handlelength=1.0)

    ax.grid(False)
    ax.tick_params(length=3)
    fig.tight_layout()

    os.makedirs(OUT_DIR, exist_ok=True)
    plt.savefig(f"{OUT_DIR}/Fig_Procrustes.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUT_DIR}/Fig_Procrustes.png", format="png", dpi=300, bbox_inches="tight")
    print(f"Procrustes figure saved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
