# Stereo-cell transcriptomic analysis of the mouse cochlea

Key analysis scripts for the mouse cochlear Stereo-cell study. Stereo-cell is
distinct from Stereo-seq. This repository contains representative analysis
steps, not an automated workflow from FASTQ files to every manuscript panel.

## Organization

| Script | Analysis | Figures |
|---|---|---|
| 01_Preprocessing_Integration.R | SoupX, QC, P1 SC/SN CCA integration, scDblFinder | 1 |
| 02_Transcriptomic_Analysis.R | Method comparisons, markers, GO-BP GSEA, hearing-loss genes | 1–5 |
| 03_SCENIC_Inference.py | GRNBoost2, motif pruning, AUCell, final target export | 2–4 |
| 04_Regulon_Analysis.R | Mean AUC, relative activity, ridge plots, optional LR overlap | 2–4 |
| 05_CellChat_Analysis.R | Separate SGN–Glia and SGN–Epi analyses and representative plots | 3–4, supplements |

Run from the repository root. Supply inputs in data/ and reference databases
in references/. Outputs go to results/. Data files are not bundled.
Final biological annotations must already exist in downstream objects.

## 1. Preprocessing and integration

Representative sample preprocessing:

    Rscript 01_Preprocessing_Integration.R preprocess data/Cell_GetExp_gene_with_background.txt Y40235F9

The input Stereo-cell expression table requires geneID, label, x, y and
MIDCount columns; label 0 indicates background. Repeat with the corresponding
input and sample ID for each sample. The example preserves nFeature_RNA > 200,
nFeature_RNA < 10000, nCount_RNA < 20000 and percent.mt < 10.

Prepare data/Y40235F9-QC.rds (SN) and data/Y40236L9-QC.rds (SC), then run:

    Rscript 01_Preprocessing_Integration.R integrate

CCA uses 3,000 features and dimensions 1–30; PCA uses 30 PCs and clustering
resolution is 0.5. The original scDblFinder source grouping is retained.
Output: results/01/P1_singlets.rds. Add final celltype annotations separately
and save as data/P1_annotated.rds.

Adult Stereo-cell/10x integration is not independently reproduced by this
P1 example; module 02 takes the final annotated adult object as input.

## 2. Transcriptomic analyses

Run independent tasks as needed:

    Rscript 02_Transcriptomic_Analysis.R p1
    Rscript 02_Transcriptomic_Analysis.R adult
    Rscript 02_Transcriptomic_Analysis.R subtype SGN
    Rscript 02_Transcriptomic_Analysis.R subtype Glia
    Rscript 02_Transcriptomic_Analysis.R subtype Epi
    Rscript 02_Transcriptomic_Analysis.R hearing P1
    Rscript 02_Transcriptomic_Analysis.R hearing Adult

Required example inputs:

| Task | File | Metadata |
|---|---|---|
| p1 | data/P1_annotated.rds | celltype; source = Single-Nuclear / Single-Cell |
| adult | data/Merge-Cochlea-Final-V3.rds | source = snRNA / StereoCell |
| subtype SGN | data/Merge-SGN-Final.rds | subtype_Final |
| subtype Glia | data/Merge-Glia-Final.rds | Glia_subtype |
| subtype Epi | data/Merge-Epi-Final-V2.rds | subtype_Final |
| hearing P1 | data/P1_SN_annotated.rds | celltype |
| hearing Adult | data/Adult_StereoCell_annotated.rds | celltype |

If final SGN identities are stored in subtype_Final_new, adjust the group
field before running modules 02 and 04.

The hearing task additionally requires data/Hearing_loss_genes_HL-adult-P1.csv
with P1, P1_Disease, Adult and Adult_Disease columns.

- P1: per-type/source counts, mean genes/UMIs, source composition and cumulative
  detected genes (nonzero in at least one cell/nucleus).
- Adult: gene overlap at the 1% detection threshold and Pearson correlation/
  linear regression of log1p-transformed mean expression.
- Subtypes: marker tables, representative UMAP, mean-expression correlation
  and GO-BP GSEA. Cell-type-specific marker thresholds from the notebooks
  are retained rather than replaced with a common threshold.
- GSEA: signed avg_log2FC ranking, logfc.threshold=0, min.pct=0.01,
  return.thresh=Inf and gene set sizes 15–500.
- Hearing genes: independent stage-specific means and row Z-scores. Symbol
  capitalization follows the supplied figure code, not an orthology lookup.

## 3. SCENIC inference

Place these resources in references/:

    mm_mgi_tfs.txt
    mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
    mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
    motifs-v9-nr.mgi-m0.001-o0.0.tbl

Run each lineage separately:

    python 03_SCENIC_Inference.py --cell-type SGN --workers 8
    python 03_SCENIC_Inference.py --cell-type Glia --workers 8
    python 03_SCENIC_Inference.py --cell-type Epi --workers 8

Default input: data/Merge-{cell_type}-Final-pySCENIC.h5ad.
Use --input PATH to override. The h5ad must contain raw counts in
layers["counts"] and unlogged counts in X for the normalization/HVG block.

The notebook strategy is retained: infer the network using count columns
selected by HVG analysis, and use the full input count matrix for module
construction and AUCell. Worker count is configurable.

Outputs in results/03/{cell_type}/:

- {cell_type}_adjacencies.csv: coexpression edges.
- motifs.csv: motif-pruning results.
- regulons.p: final regulons.
- {cell_type}_auc_matrix.csv: cells × regulons.
- {cell_type}_regulon_targets.csv: final motif-pruned target membership.

Use the final target table, not raw adjacency edges, for TF–target associations.

## 4. Regulon activity and communication links

    Rscript 04_Regulon_Analysis.R SGN
    Rscript 04_Regulon_Analysis.R Glia
    Rscript 04_Regulon_Analysis.R Epi

AUC rows are aligned to metadata by barcode. Full mean-AUC and row-Z-score
tables are exported before selecting regulons for figures. Heatmap lists follow
the supplied notebook; ridge lists follow the specified manuscript panels.
Edit the lists if final panel selections change. Missing selected regulons
trigger an error instead of being silently removed.

Optional RSS is disabled. Enable run_rss only if the R package SCENIC is
installed; its version was not supplied for this code package.

For optional TF–CellChat gene overlap, supply:

    data/SGN_Glia_bidirectional_LR.csv
    data/SGN_Epi_bidirectional_LR.csv

These may be copied from module 05 outputs or replaced by the existing final
LR tables. Set cell labels to match the tables. Ligand/receptor complex
components are matched to final regulon targets on the appropriate sender/
receiver side. SGN–Glia and SGN–Epi remain separate.
Overlap indicates candidate association, not causal or subtype-specific regulation.

## 5. CellChat

    Rscript 05_CellChat_Analysis.R SGN_Glia data/Merge-SGN_GLia.rds subtype_1004
    Rscript 05_CellChat_Analysis.R SGN_Epi data/Merge-SGN_Epi.rds subtype_1004

Combined input objects must contain normalized RNA data, SGNs and the relevant
partner cells. The third argument is the metadata grouping column. Edit the
partner labels in the script to match the objects.

The mouse database, raw.use=TRUE and min.cells=10 are retained. Sequential
execution is used for portability. raw.use refers to unprojected signaling
expression in CellChat, not raw UMI counts.

The approximately strongest 50% of positive aggregate edges are retained only
for network display. LR display tables use prob > 0 and pval < 0.05, without
the 50% aggregate-edge filter. Results are stored separately in results/05/.

## Software

Versions reported for the study:

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
| pyscenic | 0.12.1 |
| arboreto | 0.1.6 |
| ctxcore | 0.2.0 |
| scanpy | 1.9.8 |

Additional dependencies: Matrix, data.table, tidyr, future, ComplexHeatmap,
numpy, pandas, scipy and dask. This is not a fully pinned environment.

## Scope and validation

This compact release retains core computations and representative plots.
Raw-read processing, manual annotation and complete panel-by-panel figure
reconstruction are outside its scope. R and Python syntax checks passed;
the biological analyses were not rerun during packaging.
