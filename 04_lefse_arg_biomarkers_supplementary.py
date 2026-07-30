"""
LEfSe ARG Biomarker Dot Chart (Supplementary Figure)
Peshawar Wastewater Metagenomics

Genus-level LEfSe was statistically significant and forms Figure 2B of the
main text. The equivalent ARG-level LEfSe analysis (this script) did not
reach a consistent significant biomarker set across sectors and was moved to
the Supplementary Material, where it is presented alongside the drug-class
colour key produced here.

Input:
    A BioinCloud LEfSe output folder containing an "*LDA2*" results file and
    a matching ARG abundance/description file ("*abundance*" or
    "*AMR.unstratified*"). Point DATA_DIR below at that folder.

Output:
    Fig_LEfSe_LDA_Final.svg / .png -- ARGs ranked by |LDA| score, coloured by
    drug class, with the LDA=2.0 significance threshold marked.

Requirements: pandas, numpy, matplotlib; optional: adjustText (pip install adjustText)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import glob

try:
    from adjustText import adjust_text
    HAS_ADJUST = True
except ImportError:
    HAS_ADJUST = False
    print("Note: adjustText not installed; label positions will not be auto-adjusted. "
          "Install via: pip install adjustText")

plt.rcParams.update({
    "font.size": 10, "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica"],
    "axes.labelsize": 11, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 8, "figure.figsize": (8, 10), "figure.dpi": 300,
    "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})

# Point this at your local BioinCloud LEfSe output folder for the ARG-level run
DATA_DIR = "data/LEfSe_ARG"


def find_lda_file(base_path):
    lda_file = next((f for f in glob.glob(os.path.join(base_path, "*LDA2*")) if os.path.isfile(f)), None)
    if not lda_file:
        raise FileNotFoundError("LDA2 file not found in LEfSe folder.")
    return lda_file


def find_abundance_file(start_path):
    """Search the LEfSe folder, then its parent and grandparent, for the ARG
    abundance/description table (BioinCloud's folder depth varies by run)."""
    patterns = ["*abundance*", "*AMR.unstratified*"]
    for level_path in (start_path, os.path.dirname(start_path), os.path.dirname(os.path.dirname(start_path))):
        for p in patterns:
            m = glob.glob(os.path.join(level_path, p))
            if m:
                return m[0]
    return None


def load_abundance_safe(filepath):
    try:
        df = pd.read_csv(filepath, sep="\t", header=None, dtype=str)
        if "# Gene Family" in str(df.iloc[0, 0]) or "Gene Family" in str(df.iloc[0, 0]):
            df.columns = df.iloc[0]
            df = df.iloc[1:].reset_index(drop=True)
        df.columns = [str(c).replace("# ", "").strip() for c in df.columns]
        return df
    except Exception:
        df = pd.read_excel(filepath, header=0)
        df.columns = [str(c).replace("# ", "").strip() for c in df.columns]
        return df


def infer_drug_class(desc):
    if pd.isna(desc):
        return "other"
    d = str(desc).lower()
    if "vancomycin" in d or "glycopeptide" in d:
        return "glycopeptide"
    if "fluoroquinolone" in d or "quinolone" in d or "qnr" in d:
        return "fluoroquinolone"
    if "beta-lactam" in d or "betalactam" in d or "tem" in d or "shv" in d or "oxa" in d:
        return "beta-lactam"
    if "macrolide" in d or "erm" in d or "msr" in d:
        return "macrolide"
    if "aminoglycoside" in d or "aph" in d or "aac" in d:
        return "aminoglycoside"
    if "tetracycline" in d or "tet" in d:
        return "tetracycline"
    if "sulfonamide" in d or "sul" in d:
        return "sulfonamide"
    if "efflux" in d or "emr" in d or "pump" in d:
        return "efflux"
    return "other"


CLASS_COLORS = {
    "beta-lactam": "#E63946", "fluoroquinolone": "#F4A261", "macrolide": "#2A9D8F",
    "aminoglycoside": "#457B9D", "sulfonamide": "#8338EC", "tetracycline": "#FB8500",
    "efflux": "#ADB5BD", "glycopeptide": "#1D3557", "phenicol": "#E76F51",
    "trimethoprim": "#264653", "other": "#ADB5BD",
}

# Optional: label a handful of genes of interest on the plot (ARO ID -> label)
TARGET_GENES = {
    "3002801": "QnrVC4", "3002802": "QnrVC5", "3003193": "QnrVC7",
    "3002639": "APH(3'')-Ib", "3000615": "MsrA", "3000410": "Sul1",
    "3003741": "mphE", "3000186": "Tet(M)",
}


def main():
    lda_file = find_lda_file(DATA_DIR)
    abundance_file = find_abundance_file(DATA_DIR)
    if not abundance_file:
        raise FileNotFoundError("Abundance/description file not found near the LEfSe folder.")

    print(f"LDA file: {os.path.basename(lda_file)}")
    print(f"Abundance file: {os.path.basename(abundance_file)}")

    lda_raw = pd.read_csv(lda_file, sep="\t", header=None, dtype=str)
    lda_df = lda_raw.iloc[:, :2].copy()
    lda_df.columns = ["ARO_ID", "LDA_score"]
    lda_df["LDA_score"] = pd.to_numeric(lda_df["LDA_score"].str.strip().replace("-", np.nan), errors="coerce")
    lda_df = lda_df.dropna(subset=["LDA_score"])
    lda_df["ARO_ID"] = lda_df["ARO_ID"].str.strip().str.replace("ARO:", "", regex=False).str.replace("ARO_", "", regex=False)
    lda_df["LDA_abs"] = lda_df["LDA_score"].abs()

    abundance_df = load_abundance_safe(abundance_file)

    if "Description" not in abundance_df.columns:
        print("Note: 'Description' column not found; drug class will default to 'other'.")
        abundance_df["ARO_ID_clean"] = abundance_df.iloc[:, 0].str.replace("ARO:", "", regex=False)
        mapping = pd.DataFrame(index=abundance_df["ARO_ID_clean"], columns=["drug_class"]).fillna("other")
    else:
        abundance_df["ARO_ID_clean"] = (
            abundance_df.iloc[:, 0].str.replace("ARO:", "", regex=False).str.replace("ARO_", "", regex=False)
        )
        abundance_df["drug_class"] = abundance_df["Description"].apply(infer_drug_class)
        mapping = abundance_df.set_index("ARO_ID_clean")["drug_class"]

    lda_df_sorted = lda_df.sort_values("LDA_abs", ascending=False).reset_index(drop=True)
    lda_df_sorted["drug_class"] = lda_df_sorted["ARO_ID"].map(mapping).fillna("other")

    fig, ax = plt.subplots(figsize=(8, 10))
    y_pos = np.arange(len(lda_df_sorted))
    colors = [CLASS_COLORS.get(c, CLASS_COLORS["other"]) for c in lda_df_sorted["drug_class"]]

    ax.scatter(lda_df_sorted["LDA_score"], y_pos, c=colors, s=25, alpha=0.75, edgecolors="none", zorder=2)
    ax.set_xlim(0, 5.5)
    ax.axvline(x=2.0, color="#D62728", linestyle="--", linewidth=1.5, alpha=0.8, label="LDA threshold (2.0)")
    ax.axvline(x=0, color="gray", linestyle="-", linewidth=0.5, alpha=0.5)

    texts = []
    for i, row in lda_df_sorted.iterrows():
        aro = row["ARO_ID"]
        if aro in TARGET_GENES:
            t = ax.text(row["LDA_score"] + 0.05, y_pos[i], TARGET_GENES[aro],
                        fontsize=7, fontstyle="italic", va="center", ha="left", alpha=0.9)
            texts.append(t)

    if HAS_ADJUST and texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
                    expand_points=(1.2, 3.0), lim=200, only_move={"points": "y", "texts": "y"})

    ax.set_xlabel("LDA Score (Effect Size)", fontsize=11, fontweight="bold")
    ax.set_ylabel("ARGs (Ranked by LDA)", fontsize=11, fontweight="bold")
    ax.set_ylim(-0.5, len(lda_df_sorted) - 0.5)
    ax.invert_yaxis()
    ax.grid(True, axis="x", linestyle="-", alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)

    present_classes = sorted(set(lda_df_sorted["drug_class"]))
    legend_elements = [mpatches.Patch(color=CLASS_COLORS.get(c, CLASS_COLORS["other"]), label=c)
                       for c in present_classes if c in CLASS_COLORS]
    ax.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1.02, 0.5),
              ncol=1, framealpha=0.9, fontsize=7, title="Drug Class", title_fontsize=8)

    fig.subplots_adjust(left=0.15, right=0.78, top=0.96, bottom=0.04)

    os.makedirs("output", exist_ok=True)
    out_svg = "output/Fig_LEfSe_LDA_Final.svg"
    out_png = "output/Fig_LEfSe_LDA_Final.png"
    plt.savefig(out_svg, format="svg", dpi=300, bbox_inches="tight")
    plt.savefig(out_png, format="png", dpi=300, bbox_inches="tight")
    print(f"Saved: {out_svg}, {out_png}")

    print("\nDrug class distribution:")
    print(lda_df_sorted["drug_class"].value_counts())


if __name__ == "__main__":
    main()
