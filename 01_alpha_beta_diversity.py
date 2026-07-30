"""
Alpha and Beta Diversity Calculations
Peshawar Wastewater Metagenomics (hospital / slaughterhouse / community)

Reproduces the alpha diversity (Shannon, Simpson) and within/between-sector
Bray-Curtis distance values reported in the Results (bacterial community
composition and beta-diversity sections; Figure 1B).

Input:
    All_Taxa_OTU.xls   Tab-separated genus/species-level read counts per sample
    BRAY_C_1.XLS       Pre-computed Bray-Curtis dissimilarity matrix (BioinCloud output)

    PCoA ordination (axis % variance explained) and the PERMANOVA result shown in
    Figure 1B were produced directly by the BioinCloud platform (see Methods,
    Wekemo BioinCloud) and are not recomputed here.

Output:
    Shannon index, Simpson (1-D), and species richness per sample; sector means
    and ranges; intra- and inter-sector Bray-Curtis distances.
"""

import math
from collections import defaultdict

# ---------------------------------------------------------------------------
# 1. File paths (place these files in a local "data/" folder next to this script)
# ---------------------------------------------------------------------------
OTU_FILE = "data/All_Taxa_OTU.xls"
BRAY_FILE = "data/BRAY_C_1.XLS"

SAMPLES = ["PCW1", "PCW2", "PHW1", "PHW2", "PSLW1", "PSLW2"]
SECTORS = {"PHW": ["PHW1", "PHW2"], "PSLW": ["PSLW1", "PSLW2"], "PCW": ["PCW1", "PCW2"]}


# ---------------------------------------------------------------------------
# 2. Alpha diversity functions
# ---------------------------------------------------------------------------
def shannon(counts):
    """Shannon entropy H' = -sum(p_i * ln(p_i))"""
    total = sum(counts)
    if total == 0:
        return 0.0
    h = 0.0
    for c in counts:
        if c > 0:
            p = c / total
            h -= p * math.log(p)
    return h


def simpson_1_minus_d(counts):
    """Simpson diversity index, reported as 1 - D (higher = more diverse)."""
    total = sum(counts)
    if total <= 1:
        return 0.0
    d = sum(c * (c - 1) for c in counts) / (total * (total - 1))
    return 1 - d


def richness(counts):
    """Number of taxa with non-zero abundance."""
    return sum(1 for c in counts if c > 0)


def load_otu_table(path):
    """Load a tab-separated OTU table into {sample: [counts]} keyed by taxon rows."""
    data = defaultdict(list)
    with open(path, encoding="utf-8", errors="ignore") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        sample_idx = {s: header.index(s) for s in SAMPLES if s in header}
        for line in fh:
            row = line.rstrip("\n").split("\t")
            for s, idx in sample_idx.items():
                try:
                    data[s].append(float(row[idx]))
                except (IndexError, ValueError):
                    data[s].append(0.0)
    return data


def load_bray_curtis_matrix(path):
    """Load a pre-computed Bray-Curtis dissimilarity matrix (samples x samples)."""
    with open(path, encoding="utf-8", errors="ignore") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        labels = header[1:]
        matrix = {}
        for line in fh:
            row = line.rstrip("\n").split("\t")
            s1 = row[0]
            matrix[s1] = {labels[i]: float(v) for i, v in enumerate(row[1:])}
    return matrix


# ---------------------------------------------------------------------------
# 3. Run
# ---------------------------------------------------------------------------
def main():
    otu = load_otu_table(OTU_FILE)

    print("=== ALPHA DIVERSITY PER SAMPLE ===")
    alpha_results = {}
    for s in SAMPLES:
        counts = otu.get(s, [])
        h = shannon(counts)
        simp = simpson_1_minus_d(counts)
        rich = richness(counts)
        alpha_results[s] = {"shannon": h, "simpson": simp, "richness": rich}
        print(f"  {s:6s}  Shannon={h:.3f}  Simpson(1-D)={simp:.3f}  Richness={rich}")

    print("\n=== SECTOR MEANS ===")
    for sector, members in SECTORS.items():
        h_vals = [alpha_results[m]["shannon"] for m in members]
        s_vals = [alpha_results[m]["simpson"] for m in members]
        print(f"  {sector}: Shannon mean={sum(h_vals)/len(h_vals):.3f} "
              f"(range {min(h_vals):.2f}-{max(h_vals):.2f}), "
              f"Simpson mean={sum(s_vals)/len(s_vals):.3f}")

    print("\n=== BRAY-CURTIS DISTANCES (from pre-computed matrix) ===")
    bc = load_bray_curtis_matrix(BRAY_FILE)

    print("Intra-sector (within-sector, temporal replicate pairs):")
    for sector, members in SECTORS.items():
        if len(members) == 2:
            d = bc[members[0]][members[1]]
            print(f"  {sector} ({members[0]} vs {members[1]}): BC distance = {d:.3f}")

    print("\nInter-sector (all cross-sector sample pairs):")
    inter_vals = []
    for i, s1 in enumerate(SAMPLES):
        for s2 in SAMPLES[i + 1:]:
            sec1 = next(sec for sec, members in SECTORS.items() if s1 in members)
            sec2 = next(sec for sec, members in SECTORS.items() if s2 in members)
            if sec1 != sec2:
                inter_vals.append(bc[s1][s2])
    if inter_vals:
        print(f"  Mean inter-sector BC distance = {sum(inter_vals)/len(inter_vals):.3f}")


if __name__ == "__main__":
    main()
