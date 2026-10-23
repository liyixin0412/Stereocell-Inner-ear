# Stereo-cell analysis of the mouse cochlea

Core analysis scripts for mouse cochlear Stereo-cell data. The scripts retain
the study's principal analysis parameters and include basic plotting examples.

## Scripts

| File | Analysis and example plots |
|---|---|
| 01_Preprocessing_Integration.R | SoupX correction, QC, CCA integration and scDblFinder; QC violins and UMAP |
| 02_Transcriptomic_Analysis.R | P1 SC/SN and adult Stereo-cell/10x comparisons; subtype markers and GO-BP GSEA; violins, proportions, Venn, correlations, marker dots and enrichment bars |
| 03_SCENIC_Inference.py | Lineage-specific GRNBoost2, motif pruning and AUCell |
| 04_Regulon_Analysis.R | Mean AUC, relative activity, heatmaps and an example ridge plot |
| 05_CellChat_Analysis.R | CellChat inference, result tables, a network, a pathway heatmap and an example ligand–receptor bubble plot |

## Use

Set the repository directory as the working directory. Open the required script,
edit its input paths and annotation column, and run the relevant sections.
R example calls are at the end of each script; Python settings are at the top.
There is no command-line task selector or automatic whole-study driver.

- Put expression tables and annotated objects in `data/`, SCENIC databases in
  `references/`; outputs are written to `results/`.
- **01:** repeat sample preprocessing for each Stereo-cell sample. Run integration
  only after all required QC objects, including the public adult 10x object, exist.
- **02:** method comparisons use final annotated P1/adult objects. The subtype
  section uses Glia as an example; repeat with the corresponding SGN/Epi object.
- **03 → 04:** run separately for SGN, Glia and Epi. The h5ad must contain raw counts
  in `layers["counts"]`; `X` must be suitable for normalization and log transformation.
  Match AUC cell barcodes to the annotated Seurat object. For 03, activate the
  study Python environment and run `python 03_SCENIC_Inference.py`.
- **05:** use SGN–Glia and SGN–Epi as separate input objects. An existing CellChat
  object can be used for plotting without rerunning inference.

Example filenames refer to the study's processed objects. Data, manual biological
annotations and reference databases are supplied separately; update paths to
match local files. These scripts do not start from FASTQ files.

SCENIC reference files:

    mm_mgi_tfs.txt
    mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
    mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
    motifs-v9-nr.mgi-m0.001-o0.0.tbl

## Software

The following software versions were used in the study.

| Software | Version |
|---|---|
| Seurat | 4.3.0 |
| SeuratObject | 4.1.3 |
| SoupX | 1.6.2 |
| scDblFinder | 1.12.0 |
| SingleCellExperiment | 1.20.1 |
| msigdbr | 7.5.1 |
| fgsea | 1.24.0 |
| CellChat | 1.6.1 |
| pheatmap | 1.0.12 |
| ggridges | 0.5.6 |
| ggplot2 | 3.5.1 |
| dplyr | 1.1.4 |
| Python | 3.8.13 |
| pySCENIC | 0.12.1 |
| arboreto | 0.1.6 |
| ctxcore | 0.2.0 |
| scanpy | 1.9.8 |

Python and the SCENIC packages are used in module 03. The remaining analysis
and visualization modules use R.

Additional dependencies (versions not specified here):

- Matrix
- data.table
- tidyr
- patchwork
- ggrepel
- VennDiagram
- viridisLite
- future
- ComplexHeatmap
- numpy
- pandas
- scipy
- dask
- distributed

Use the original analysis environment where possible. The version list documents
the study software; it is not a complete environment lockfile.

## Scope

This repository provides the core computational workflow and illustrative plots.
It does not include the hearing-loss gene analysis, manual annotation decisions,
manuscript-specific marker/TF/pathway selections or final figure assembly.
Example markers and GSEA terms are ranked from the results; four variable regulons
and one inferred CellChat pathway illustrate plotting. These display choices do
not constitute biological-priority selection or reproduce the final figure panels.

The scripts were syntax-checked and selected calculations and plots were tested
with synthetic data. Full inference was not rerun on the study data during code
organization; numerical identity with the published figures is not guaranteed.
