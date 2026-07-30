"""
Integrated Resistome-Mobilome Spearman Correlation Analysis
Peshawar Wastewater Metagenomics -- Figure 6 (Integrated correlation heatmaps)

Computes Spearman correlations between (1) bacterial genus abundance and ARG
gene family abundance, and (2) MGE type abundance and ARG gene family
abundance, across the six samples. Reproduces the rho/p-value matrices behind
Figure 6A (taxa-ARG) and Figure 6B (MGE-ARG).

Input:
    OTU_TA_3.XLS               Genus-level relative abundance, tab-separated
    All_CARD_AMR_Gene_Family.txt  ARG gene family abundance, tab-separated
    All_MGE_Level1.txt          MGE Level-1 type abundance, tab-separated

Output:
    output/rho_taxa_arg.csv   15 genera x 12 ARG families Spearman rho matrix
    output/rho_mge_arg.csv    10 MGE types x 12 ARG families Spearman rho matrix
"""

import os
import pandas as pd
import numpy as np
from scipy import stats

SAMPLES = ["PCW1", "PCW2", "PHW1", "PHW2", "PSLW1", "PSLW2"]
DATA_DIR = "data"
OUT_DIR = "output"


def load_data():
    tax = pd.read_csv(f"{DATA_DIR}/OTU_TA_3.XLS", sep="\t")
    tax = tax[["name"] + SAMPLES]
    tax["total"] = tax[SAMPLES].sum(axis=1)

    arg = pd.read_csv(f"{DATA_DIR}/All_CARD_AMR_Gene_Family.txt", sep="\t")
    arg = arg.rename(columns={"# CARD AMR Gene Family": "Gene_Family"})
    arg = arg[["Gene_Family"] + SAMPLES]
    arg["total"] = arg[SAMPLES].sum(axis=1)

    mge = pd.read_csv(f"{DATA_DIR}/All_MGE_Level1.txt", sep="\t")
    mge = mge.rename(columns={"# MGE Level1": "MGE"})
    mge = mge[["MGE"] + SAMPLES]
    mge["total"] = mge[SAMPLES].sum(axis=1)

    return tax, arg, mge


def spearman_matrix(rows_df, row_label_col, cols_df, col_label_col):
    """Pairwise Spearman rho/p between every row feature and every column feature."""
    row_arr = rows_df[SAMPLES].values
    col_arr = cols_df[SAMPLES].values
    rho = np.zeros((len(rows_df), len(cols_df)))
    pval = np.zeros((len(rows_df), len(cols_df)))
    for i in range(len(rows_df)):
        for j in range(len(cols_df)):
            r, p = stats.spearmanr(row_arr[i], col_arr[j])
            rho[i, j] = r
            pval[i, j] = p
    return rho, pval


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    tax, arg, mge = load_data()

    print(f"Taxonomy: {len(tax)} genera | ARG families: {len(arg)} | MGE types: {len(mge)}")

    # Select top genera (excluding Vibrio -- single-timepoint bloom artefact -- and unclassified)
    exclude_genera = {"Vibrio", "unclassified"}
    top_genera_all = tax.nlargest(30, "total")
    top_genera = top_genera_all[~top_genera_all["name"].isin(exclude_genera)].head(15)

    top_arg = arg.nlargest(12, "total")

    mge_excl = {"unclassified"}
    top_mge = mge[~mge["MGE"].isin(mge_excl)].nlargest(10, "total")

    print("\nSelected genera:", top_genera["name"].tolist())
    print("Selected ARG families:", top_arg["Gene_Family"].tolist())
    print("Selected MGE types:", top_mge["MGE"].tolist())

    # ── Taxa-ARG correlation ──────────────────────────────────────────────────
    rho_taxa_arg, p_taxa_arg = spearman_matrix(top_genera, "name", top_arg, "Gene_Family")

    pairs = []
    for i, gen in enumerate(top_genera["name"].tolist()):
        for j, fam in enumerate(top_arg["Gene_Family"].tolist()):
            pairs.append((gen, fam, rho_taxa_arg[i, j], p_taxa_arg[i, j]))
    pairs_sorted = sorted(pairs, key=lambda x: abs(x[2]), reverse=True)

    print("\n=== TOP 20 TAXA-ARG PAIRS BY |RHO| ===")
    for gen, fam, rho, p in pairs_sorted[:20]:
        print(f"  {gen:30s} x {fam:35s} rho={rho:+.3f} p={p:.3f}")

    # ── MGE-ARG correlation ───────────────────────────────────────────────────
    rho_mge_arg, p_mge_arg = spearman_matrix(top_mge, "MGE", top_arg, "Gene_Family")

    pairs_mge = []
    for i, mg in enumerate(top_mge["MGE"].tolist()):
        for j, fam in enumerate(top_arg["Gene_Family"].tolist()):
            pairs_mge.append((mg, fam, rho_mge_arg[i, j], p_mge_arg[i, j]))
    pairs_mge_sorted = sorted(pairs_mge, key=lambda x: abs(x[2]), reverse=True)

    print("\n=== TOP 20 MGE-ARG PAIRS BY |RHO| ===")
    for mg, fam, rho, p in pairs_mge_sorted[:20]:
        print(f"  {mg:30s} x {fam:35s} rho={rho:+.3f} p={p:.3f}")

    print("\n=== SUMMARY ===")
    print(f"Taxa-ARG: {np.sum(np.abs(rho_taxa_arg) > 0.8)} pairs with |rho|>0.8")
    print(f"Taxa-ARG: {np.sum(np.abs(rho_taxa_arg) > 0.9)} pairs with |rho|>0.9")
    print(f"MGE-ARG:  {np.sum(np.abs(rho_mge_arg) > 0.8)} pairs with |rho|>0.8")
    print(f"MGE-ARG:  {np.sum(np.abs(rho_mge_arg) > 0.9)} pairs with |rho|>0.9")

    pd.DataFrame(rho_taxa_arg, index=top_genera["name"].tolist(),
                 columns=top_arg["Gene_Family"].tolist()).to_csv(f"{OUT_DIR}/rho_taxa_arg.csv")
    pd.DataFrame(rho_mge_arg, index=top_mge["MGE"].tolist(),
                 columns=top_arg["Gene_Family"].tolist()).to_csv(f"{OUT_DIR}/rho_mge_arg.csv")
    print(f"\nMatrices saved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
