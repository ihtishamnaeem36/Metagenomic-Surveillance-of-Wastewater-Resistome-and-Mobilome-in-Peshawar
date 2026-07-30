"""
Dual-Carriage Genus Analysis and Bubble Plot
Peshawar Wastewater Metagenomics -- Figure 9A (dual-carriage composite)

Identifies genera that carry both genus-stratified ARG signal (CARD) and
genus-stratified MGE signal, and plots the top ARG/MGE dual carriers as a
bubble scatter (x = MGE total, y = ARG total; bubble size = number of ARG
features; colour = dominant ARG sector; edge thickness = Spearman rho with
the paired MGE profile).

Input:
    CARD_abundance_mixed.xls   Genus-stratified CARD (ARG) abundance table
    MGE_abundance_mixed.xls    Genus-stratified MGE abundance table

    Both files use the format "ARO_or_MGE_ID|k__...;g__Genus.s__..." for
    genus-stratified rows.

Intermediate step not automated here: per-genus ARG_total, MGE_total,
Spearman rho, and dominant ARG sector (used in the plotting step below) were
computed by aggregating each dual-carrying genus's per-sample ARG and MGE
totals and correlating them across the six samples (the same approach as
scripts/python/06_resistome_mobilome_spearman_correlation.py, applied per
genus rather than per feature). Save that intermediate table as
data/dual_carrying_genera.csv with columns:
genus, ARG_total, MGE_total, rho, p_rho, dom_ARG_sector, n_ARG_features

Output:
    Fig_DualCarriage.svg / .png (300 dpi)
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = "data"
OUT_DIR = "output"
SAMPLES = ["PCW1", "PCW2", "PHW1", "PHW2", "PSLW1", "PSLW2"]


def extract_genus(feature_id):
    """Extract genus from a 'ID|k__...;g__Genus.s__...' stratified feature string."""
    m = re.search(r"g__([^.]+)\.s__", feature_id)
    if m:
        g = m.group(1).strip()
        return g if g else None
    return None


def find_dual_carrying_genera():
    """Identify genera with both genus-stratified ARG and MGE signal."""
    card = pd.read_csv(f"{DATA_DIR}/CARD_abundance_mixed.xls", sep="\t")
    card = card.rename(columns={"CARD": "feature"})
    card_strat = card[card["feature"].str.contains(r"\|", na=False)].copy()
    card_strat = card_strat[~card_strat["feature"].str.contains("unclassified", case=False, na=False)]
    card_strat["genus"] = card_strat["feature"].apply(extract_genus)
    card_strat = card_strat[card_strat["genus"].notna() & (card_strat["genus"] != "")]

    mge = pd.read_csv(f"{DATA_DIR}/MGE_abundance_mixed.xls", sep="\t")
    mge = mge.rename(columns={"MGE": "feature"})
    mge_strat = mge[mge["feature"].str.contains(r"\|", na=False)].copy()
    mge_strat = mge_strat[~mge_strat["feature"].str.contains("unclassified", case=False, na=False)]
    mge_strat["genus"] = mge_strat["feature"].apply(extract_genus)
    mge_strat = mge_strat[mge_strat["genus"].notna() & (mge_strat["genus"] != "")]

    print(f"CARD stratified rows: {len(card_strat)}  |  Unique genera: {card_strat['genus'].nunique()}")
    print(f"MGE stratified rows:  {len(mge_strat)}  |  Unique genera: {mge_strat['genus'].nunique()}")

    dual = sorted(set(card_strat["genus"].unique()) & set(mge_strat["genus"].unique()))
    print(f"\nDual-carrying genera (ARG + MGE): {len(dual)}")
    print(dual)
    return card_strat, mge_strat, dual


def plot_dual_carriage_bubble():
    """Bubble plot of the top dual-carrying genera. Requires
    data/dual_carrying_genera.csv -- see module docstring."""
    df = pd.read_csv(f"{DATA_DIR}/dual_carrying_genera.csv")

    # Vibrio and Salmonella excluded: single-timepoint bloom artefact (see Methods)
    exclude = {"Vibrio", "Salmonella"}
    df_plot = df[~df["genus"].isin(exclude)].nlargest(15, "ARG_total").copy()
    print(f"Plotting {len(df_plot)} genera (Vibrio/Salmonella excluded)")
    print(df_plot[["genus", "ARG_total", "MGE_total", "rho", "p_rho", "dom_ARG_sector"]].to_string())

    sector_colors = {"PCW": "#1F77B4", "PHW": "#D62728", "PSLW": "#FF7F0E"}
    sector_labels = {"PCW": "Community", "PHW": "Hospital", "PSLW": "Slaughterhouse"}

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 11, "axes.labelweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 300, "savefig.dpi": 300,
    })

    fig, ax = plt.subplots(figsize=(8, 7))

    for _, row in df_plot.iterrows():
        col = sector_colors.get(row["dom_ARG_sector"], "#888888")
        size = max(60, row["n_ARG_features"] * 4)
        rho = row["rho"]
        lw = 1.0 + 2.5 * max(0, rho)  # edge thickness encodes coupling strength

        ax.scatter(row["MGE_total"], row["ARG_total"], s=size, color=col, alpha=0.82,
                   edgecolors="white" if rho < 0.8 else "#222222", linewidths=lw, zorder=4)
        ax.annotate(f"$\\it{{{row['genus']}}}$", xy=(row["MGE_total"], row["ARG_total"]),
                    xytext=(6, 4), textcoords="offset points", fontsize=8, color="#222222", zorder=5)

    ax.set_xlabel("Total MGE Abundance (CPC, stratified)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Total ARG Abundance (CPC, stratified)", fontsize=11, fontweight="bold")
    ax.set_xscale("log")
    ax.set_yscale("log")

    for n_feat, label in [(5, "5 ARG features"), (20, "20 ARG features"), (50, "50 ARG features")]:
        ax.scatter([], [], s=max(60, n_feat * 4), color="#AAAAAA", alpha=0.8,
                   edgecolors="white", label=label)
    for sec, col in sector_colors.items():
        ax.scatter([], [], s=80, color=col, alpha=0.85, edgecolors="white",
                   label=f"Dominant: {sector_labels[sec]}")
    ax.scatter([], [], s=80, color="#AAAAAA", edgecolors="#222222", linewidths=3.5,
               label="Strong coupling (rho >= 0.8)")

    ax.legend(loc="upper left", fontsize=8, framealpha=0.9, edgecolor="#CCCCCC",
              title="Legend", title_fontsize=8, ncol=1, handlelength=1.2)
    ax.grid(True, which="both", linestyle=":", alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)

    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    plt.savefig(f"{OUT_DIR}/Fig_DualCarriage.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUT_DIR}/Fig_DualCarriage.png", format="png", dpi=300, bbox_inches="tight")
    print(f"\nFigure saved to {OUT_DIR}/")


if __name__ == "__main__":
    find_dual_carrying_genera()
    plot_dual_carriage_bubble()
