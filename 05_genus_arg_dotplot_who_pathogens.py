#!/usr/bin/env python3
"""
Genus-ARG Summary Analysis
Peshawar Wastewater Metagenomics -- Figure 4A (WHO priority pathogen dot plot)

Summarises genus-level ARG abundance (mean TPM per sector) that underlies the
WHO priority pathogen-resistance gene co-occurrence dot plot, including
Kruskal-Wallis significance testing and sector-enrichment breakdowns.

Input (produced upstream in R -- see scripts/r/genus_arg_kruskal_lefse_heatmap.R):
    Genus_ARG_mean_matrix_all.csv      Feature (Genus|ARO) x sector mean TPM matrix
    Genus_ARG_mean_TPM_per_sector.csv  Long-format mean TPM per sector
    kruskal_pvalues_fdr.csv            Genus-ARO pair, raw p, FDR-adjusted p
"""

import pandas as pd
import numpy as np

DATA_DIR = "data"


def main():
    mat = pd.read_csv(f"{DATA_DIR}/Genus_ARG_mean_matrix_all.csv")
    kw = pd.read_csv(f"{DATA_DIR}/kruskal_pvalues_fdr.csv")

    print("=== MATRIX SHAPE ===")
    print(mat.shape)

    print("\n=== TOTAL TPM PER SECTOR ===")
    total_tpm = mat[["Community", "Hospital", "Slaughterhouse"]].sum()
    print(total_tpm)
    print(f"\nGrand total TPM: {total_tpm.sum():.1f}")

    # Top 20 features by maximum mean TPM
    mat["max"] = mat[["Community", "Hospital", "Slaughterhouse"]].max(axis=1)
    mat["total"] = mat[["Community", "Hospital", "Slaughterhouse"]].sum(axis=1)
    top20 = mat.nlargest(20, "max")
    print("\n=== TOP 20 BY MAX MEAN TPM ===")
    print(top20[["Feature", "Community", "Hospital", "Slaughterhouse", "max"]].to_string())

    # Genus totals (summing all AROs per genus)
    mat["Genus"] = mat["Feature"].str.split("|").str[0]
    genus_totals = mat.groupby("Genus")[["Community", "Hospital", "Slaughterhouse"]].sum()
    genus_totals["Total"] = genus_totals.sum(axis=1)
    genus_totals = genus_totals.sort_values("Total", ascending=False)
    print("\n=== GENUS TOTALS ===")
    print(genus_totals.head(20).to_string())

    # Kruskal-Wallis statistics
    print("\n=== KRUSKAL-WALLIS STATS ===")
    print(f"Total pairs tested: {len(kw)}")
    print(f"Min raw p_value: {kw['p_value'].min():.4f}")
    print(f"Min adjusted p (FDR): {kw['p_adjusted'].min():.4f}")
    print(f"Significant after FDR (adj p < 0.05): {(kw['p_adjusted'] < 0.05).sum()}")
    print("\nTop 10 lowest adjusted p:")
    print(kw.nsmallest(10, "p_adjusted")[["Genus", "ARO", "p_value", "p_adjusted"]].to_string())

    # Clinically relevant genera (WHO priority pathogens + common wastewater genera)
    clinical = ["Escherichia", "Klebsiella", "Acinetobacter", "Staphylococcus",
                "Enterococcus", "Pseudomonas", "Salmonella", "Vibrio",
                "Aliarcobacter", "Streptococcus", "Clostridium", "Bifidobacterium"]
    print("\n=== CLINICALLY RELEVANT GENERA ===")
    for g in clinical:
        subset = mat[mat["Genus"] == g]
        if len(subset) > 0:
            c = subset["Community"].sum()
            h = subset["Hospital"].sum()
            sl = subset["Slaughterhouse"].sum()
            print(f"{g}: Community={c:.2f}, Hospital={h:.2f}, Slaughterhouse={sl:.2f}, n_ARGs={len(subset)}")

    # Sector-enriched genera
    hosp_enrich = genus_totals[
        (genus_totals["Hospital"] > genus_totals["Community"]) &
        (genus_totals["Hospital"] > genus_totals["Slaughterhouse"])
    ].copy()
    hosp_enrich["H_to_C_ratio"] = hosp_enrich["Hospital"] / (hosp_enrich["Community"] + 0.001)
    print("\n=== HOSPITAL-ENRICHED GENERA (H > C and H > Slaughterhouse) ===")
    print(hosp_enrich.sort_values("Hospital", ascending=False).head(15).to_string())

    sl_enrich = genus_totals[
        (genus_totals["Slaughterhouse"] > genus_totals["Community"]) &
        (genus_totals["Slaughterhouse"] > genus_totals["Hospital"])
    ]
    print("\n=== SLAUGHTERHOUSE-ENRICHED GENERA ===")
    print(sl_enrich.sort_values("Slaughterhouse", ascending=False).head(15).to_string())

    # Specific high-priority ARGs: NDM, VIM, OXA-23, vanR
    priority_aros = ["ARO:3000027", "ARO:3000619", "ARO:3000780", "ARO:3003923"]
    print("\n=== SPECIFIC HIGH-PRIORITY ARO IDs ===")
    for aro in priority_aros:
        sub = mat[mat["Feature"].str.contains(aro, na=False)]
        if len(sub) > 0:
            print(sub[["Feature", "Community", "Hospital", "Slaughterhouse"]].to_string())

    print("\nAll analyses completed.")


if __name__ == "__main__":
    main()
