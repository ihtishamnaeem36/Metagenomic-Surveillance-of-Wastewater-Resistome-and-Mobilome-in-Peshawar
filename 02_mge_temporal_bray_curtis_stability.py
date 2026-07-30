"""
MGE Temporal Stability (Bray-Curtis Similarity)
Peshawar Wastewater Metagenomics

Computes the within-sector Bray-Curtis similarity between the two sampling
timepoints (T1, T2) for MGE Level-2 family profiles. Reproduces the temporal
stability values reported in the Results (MGE temporal dynamics section) and
underlying Figure 5B/C (MGE PCoA and slope plot panels).

Input:
    All_MGE_Level2.txt   Tab-separated MGE family abundance table (TPM), one
                         row per MGE family, one column per sample.

Output:
    Bray-Curtis similarity (1 - BC distance) for each sector's T1-vs-T2 pair.
"""

import pandas as pd
import numpy as np

SAMPLES = ["PCW1", "PCW2", "PHW1", "PHW2", "PSLW1", "PSLW2"]


def bray_curtis_similarity(a, b):
    """Bray-Curtis similarity = 1 - Bray-Curtis distance.

    BC distance = sum(|x_i - y_i|) / sum(x_i + y_i)
    """
    a, b = np.array(a), np.array(b)
    bc_dist = np.sum(np.abs(a - b)) / np.sum(a + b)
    return 1 - bc_dist


def main():
    mge = pd.read_csv("data/All_MGE_Level2.txt", sep="\t")
    mge = mge.rename(columns={"# MGE Level2": "MGE"})
    mge = mge[["MGE"] + SAMPLES]

    pairs = [
        ("PHW1", "PHW2", "Hospital"),
        ("PSLW1", "PSLW2", "Slaughterhouse"),
        ("PCW1", "PCW2", "Community"),
    ]

    print("=== BRAY-CURTIS SIMILARITY (MGE Level 2, within-sector temporal pairs) ===")
    for s1, s2, sector in pairs:
        sim = bray_curtis_similarity(mge[s1].values, mge[s2].values)
        print(f"  {sector} ({s1} vs {s2}): BC similarity = {sim:.4f}")


if __name__ == "__main__":
    main()
