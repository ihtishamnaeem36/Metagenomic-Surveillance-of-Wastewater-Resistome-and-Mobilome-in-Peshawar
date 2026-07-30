"""
LEfSe LDA Score Bar Chart
Peshawar Wastewater Metagenomics — Figure 2B

Top genus-level biomarkers per sector (LDA > 2.0), plotted as a grouped
horizontal bar chart. Reproduces Figure 2B.

Input:
    Hard-coded LDA scores below (top 15 genera per sector), extracted from the
    LEfSe output produced by the BioinCloud platform (see Methods).

Output:
    Fig2B_LEfSe_LDA_BarChart.svg / .pdf / .png (300 dpi)

Requirements: matplotlib, numpy
    pip install matplotlib numpy
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Data: top 15 genera per sector, LDA scores ────────────────────────────────
# Format: (genus, LDA_score)
hospital_data = [
    ("Thauera", 5.448), ("Aliarcobacter", 4.964), ("Diaphorobacter", 4.766),
    ("Acidovorax", 4.593), ("Quatrionicoccus", 4.510), ("Desulfovibrio", 4.380),
    ("Alicycliphilus", 4.340), ("Bifidobacterium", 4.275), ("Micropruina", 4.097),
    ("Hydrogenophaga", 4.033), ("Pseudothauera", 3.844), ("Ferribacterium", 3.814),
    ("Tessaracoccus", 3.735), ("Desulfobulbus", 3.732), ("Dechloromonas", 3.725),
]

slaughterhouse_data = [
    ("Flavobacterium", 5.004), ("Acinetobacter", 4.942), ("Acetoanaerobium", 4.888),
    ("Cloacibacterium", 4.828), ("Aeromonas", 4.565), ("Sphaerotilus", 4.562),
    ("Arcobacter", 4.482), ("Bacteroides", 4.452), ("Comamonas", 4.382),
    ("Enterococcus", 4.348), ("Shewanella", 4.235), ("Chryseobacterium", 4.169),
    ("Streptococcus", 4.167), ("Escherichia", 4.163), ("Fusobacterium", 4.110),
]

community_data = [
    ("Vibrio", 5.566), ("Pseudoxanthomonas", 4.710), ("Pseudomonas", 4.700),
    ("Mycolicibacterium", 4.303), ("Thermomonas", 4.224), ("Malaciobacter", 4.222),
    ("Stenotrophomonas", 4.032), ("Paracoccus", 4.011), ("Corynebacterium", 3.996),
    ("Phnomibacter", 3.987), ("Afipia", 3.749), ("Runella", 3.713),
    ("Rhodobacter", 3.696), ("Draconibacterium", 3.656), ("Fuscovulum", 3.569),
]

# ── Colours (colour-blind friendly) ───────────────────────────────────────────
COLORS = {"Hospital": "#D62728", "Slaughterhouse": "#FF7F0E", "Community": "#1F77B4"}
EDGE_COLORS = {"Hospital": "#8B0000", "Slaughterhouse": "#B85800", "Community": "#0A4A7A"}

sections = [
    ("Hospital", hospital_data),
    ("Slaughterhouse", slaughterhouse_data),
    ("Community", community_data),
]


def main():
    fig_height, fig_width = 14, 9
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    y_pos, y_labels, y_colors, y_lda, y_edge = [], [], [], [], []
    current_y = 0
    sector_label_y, sector_spans = {}, {}

    for sector, data in sections:
        start_y = current_y
        for genus, lda in reversed(data):  # highest LDA plotted at top of sector
            y_pos.append(current_y)
            y_labels.append(genus)
            y_colors.append(COLORS[sector])
            y_lda.append(lda)
            y_edge.append(EDGE_COLORS[sector])
            current_y += 1
        sector_label_y[sector] = (start_y + current_y - 1) / 2
        sector_spans[sector] = (start_y - 0.5, current_y - 0.5)
        current_y += 2  # gap between sectors

    bar_height = 0.72
    for y, lda, color, edge in zip(y_pos, y_lda, y_colors, y_edge):
        ax.barh(y, lda, height=bar_height, color=color, edgecolor=edge,
                linewidth=0.5, alpha=0.88, zorder=3)
        ax.text(lda + 0.04, y, f"{lda:.2f}", va="center", ha="left",
                fontsize=7.5, color="#333333", zorder=4)

    for sector, (y0, y1) in sector_spans.items():
        ax.axhspan(y0, y1, color=COLORS[sector], alpha=0.04, zorder=1)

    x_max = max(y_lda) + 0.8
    for sector, cy in sector_label_y.items():
        ax.text(x_max + 0.15, cy, sector, va="center", ha="left",
                fontsize=10, fontweight="bold", color=COLORS[sector],
                rotation=90, rotation_mode="anchor")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=8.5, fontstyle="italic")
    ax.set_xlabel("LDA Score (log10)", fontsize=11, labelpad=8)
    ax.set_xlim(0, x_max + 1.2)
    ax.xaxis.set_tick_params(labelsize=9)
    ax.axvline(x=2.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6,
               zorder=2, label="LDA threshold = 2.0")
    ax.xaxis.grid(True, linestyle=":", linewidth=0.5, alpha=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title("LEfSe -- Top 15 Genus-Level Biomarkers per Sector\n"
                 "(LDA Score > 2.0, Bacteria-only)", fontsize=12, fontweight="bold", pad=12)

    legend_patches = [
        mpatches.Patch(facecolor=COLORS[s], edgecolor=EDGE_COLORS[s],
                       linewidth=0.8, label=s, alpha=0.88)
        for s in ["Hospital", "Slaughterhouse", "Community"]
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=9,
              framealpha=0.9, edgecolor="#CCCCCC", handlelength=1.2, handleheight=0.9)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.invert_yaxis()

    plt.tight_layout(pad=1.5)

    fig.savefig("output/Fig2B_LEfSe_LDA_BarChart.svg", format="svg", bbox_inches="tight", dpi=300)
    fig.savefig("output/Fig2B_LEfSe_LDA_BarChart.pdf", format="pdf", bbox_inches="tight", dpi=300)
    fig.savefig("output/Fig2B_LEfSe_LDA_BarChart_300dpi.png", format="png", bbox_inches="tight", dpi=300)
    print("Saved SVG, PDF, and PNG to output/")
    plt.close()


if __name__ == "__main__":
    main()
