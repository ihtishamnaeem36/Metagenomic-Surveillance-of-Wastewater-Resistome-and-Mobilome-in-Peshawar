# Peshawar Wastewater Metagenomics — Analysis Code

Code supporting the manuscript **"Metagenomic Surveillance of the Resistome and
Mobilome Across Hospital, Slaughterhouse, and Community Wastewater in Peshawar,
Pakistan"** (submitted to *Chemosphere*).

Six shotgun metagenomic wastewater samples (two hospital, two slaughterhouse,
two community, sampled 30 days apart) were processed on the Wekemo BioinCloud
platform (Gao et al., 2024) for read QC, taxonomic profiling (Kraken2/Bracken/
DIAMOND), CARD-based ARG annotation, and MGE annotation. This repository
contains the **downstream statistical analysis and figure-generation code**
written on top of BioinCloud's output tables — it is not a replacement for the
primary pipeline, which is described in the manuscript's Methods section.

## Repository structure

```
scripts/
  python/
    01_alpha_beta_diversity.py                    Shannon/Simpson diversity, Bray-Curtis distances (Fig. 1B)
    02_mge_temporal_bray_curtis_stability.py       Within-sector temporal BC similarity (Fig. 5B,C)
    03_lefse_genus_biomarkers_barplot_fig2B.py     Genus-level LEfSe bar chart (Fig. 2B)
    04_lefse_arg_biomarkers_supplementary.py       ARG-level LEfSe dot chart (Supplementary Fig. S2)
    05_genus_arg_dotplot_who_pathogens.py          Genus-ARG summary stats (Fig. 4A)
    06_resistome_mobilome_spearman_correlation.py  Taxa-ARG and MGE-ARG correlations (Fig. 6A,B)
    07_procrustes_analysis_plot.py                 Procrustes superimposition plot (Fig. 7)
    08_rda_analysis_arg_ordination.py              RDA/PERMANOVA/envfit stats + ARG biplot (Fig. 8A)
    09_dual_carriage_genus_analysis.py             Dual ARG+MGE carriage bubble plot (Fig. 9A)
  r/
    genus_arg_kruskal_lefse_heatmap.R              Genus-ARG Kruskal-Wallis/FDR + LEfSe + heatmap (Fig. 4A)
reference_outputs/
  lefse_arg_drug_class_distribution_sample_output.txt   Console output from script 04, kept for reference
requirements.txt
LICENSE
```

Figure numbers above refer to the current (post-merge) main-text numbering:
Figures 4, 5, and 9 each combine what were originally separate panels into a
single multi-panel figure. Figures 3 (circos), 5's circos panel, and parts of
Figure 6/8 were rendered directly by BioinCloud and are not reproduced by
custom code here.

## Requirements

Python packages (see `requirements.txt`):
```
pip install -r requirements.txt
```

R packages:
```r
install.packages(c("tidyverse", "reshape2", "broom", "pheatmap"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("lefser")
```

## Data

Raw sequencing data and the intermediate abundance tables each script expects
(CARD/MGE abundance tables, RDA/DCA/Procrustes output, etc.) are available from
the corresponding authors upon reasonable request (see the manuscript's Data
Availability statement). Each script's docstring lists the specific input
files it expects, normally placed in a local `data/` folder next to the
script.

## Citation

If you use this code, please cite the manuscript (citation to be added once
published) and, if relevant, this repository.

## License

MIT License — see `LICENSE`.

## Contact

Corresponding authors: Ishaq Khan (ishaq@uswat.edu.pk), Arshad Iqbal
(arshad.iqbal@uswat.edu.pk), Muhammad Shafiq (drshafiq@stu.edu.cn).
