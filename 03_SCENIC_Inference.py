"""SCENIC inference for independently prepared SGN, Glia and Epi datasets.

The supplied notebook uses layers['counts'] for network inference/AUCell
and adata.X for gene filtering, normalization and HVG selection.
"""
import glob
import pickle
from pathlib import Path


# Settings: edit for each independently prepared lineage.
cell_type = "Glia"  # SGN, Glia or Epi
input_file = f"data/Merge-{cell_type}-Final-pySCENIC.h5ad"
workers = 8
threads = 1
score_only = False  # True reuses results/03/TYPE/regulons.p


def main():
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse
    from dask.distributed import Client
    from dask.diagnostics import ProgressBar
    from arboreto.algo import grnboost2
    from arboreto.utils import load_tf_names
    from ctxcore.rnkdb import FeatherRankingDatabase
    from pyscenic.utils import modules_from_adjacencies
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell

    # 1. Input matrix and reference databases ----------------------------------
    if cell_type not in {"SGN", "Glia", "Epi"} or workers < 1 or threads < 1:
        raise ValueError("Check cell_type, workers and threads.")
    out_dir = Path("results/03") / cell_type
    out_dir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(input_file)
    if "counts" not in adata.layers:
        raise ValueError("Input h5ad must include layers['counts'].")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Cell barcodes and gene names must be unique.")
    counts = adata.layers["counts"]
    ex_matrix = pd.DataFrame(
        counts.toarray() if sparse.issparse(counts) else counts,
        index=adata.obs_names, columns=adata.var_names,
    )
    if not np.isfinite(ex_matrix.to_numpy()).all() or (ex_matrix.to_numpy() < 0).any():
        raise ValueError("Counts contain nonfinite or negative values.")
    regulon_file = out_dir / "regulons.p"

    if score_only:
        # Only load trusted local pickle files.
        with regulon_file.open("rb") as handle:
            regulons = pickle.load(handle)
    else:
        database_glob = "references/mm10__refseq-r80__*tss.mc9nr.genes_vs_motifs.rankings.feather"
        database_files = sorted(glob.glob(database_glob))
        motif_file = Path("references/motifs-v9-nr.mgi-m0.001-o0.0.tbl")
        tf_file = Path("references/mm_mgi_tfs.txt")
        if not database_files or not motif_file.exists() or not tf_file.exists():
            raise FileNotFoundError("Prepare the mouse references listed in README.")
        tf_names = load_tf_names(str(tf_file))

        # 2. HVG selection: preserve the notebook's input-X preprocessing.
        # X must be suitable for normalize_total/log1p, not already log transformed.
        sc.pp.filter_genes(adata, min_cells=10)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(
            adata, n_top_genes=5000, flavor="seurat", batch_key=None
        )
        hvg = adata.var_names[adata.var.highly_variable]
        ex_hvg = ex_matrix.loc[:, ex_matrix.columns.intersection(hvg)]
        if not set(tf_names).intersection(ex_hvg.columns):
            raise ValueError("No TFs overlap the selected HVGs; check gene symbols.")
        pd.Series(ex_hvg.columns, name="gene").to_csv(out_dir / "HVGs.csv", index=False)

        # 3. Coexpression network and motif-pruned regulons --------------------
        np.random.seed(42)
        with Client(n_workers=workers, threads_per_worker=threads) as client:
            adjacencies = grnboost2(
                expression_data=ex_hvg, tf_names=tf_names,
                client_or_address=client, verbose=True,
            )
        adjacencies.to_csv(out_dir / f"{cell_type}_adjacencies.csv", index=False)
        modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
        if not modules:
            raise ValueError("No coexpression modules retained.")
        # Original notebook exported module targets before motif pruning.
        module_targets = [
            {"Module": m.name, "TF": m.transcription_factor,
             "Target": gene, "Weight": weight}
            for m in modules for gene, weight in m.gene2weight.items()
        ]
        pd.DataFrame(module_targets).to_csv(out_dir / "unpruned_module_targets.csv", index=False)
        databases = [
            FeatherRankingDatabase(fname=f, name=Path(f).stem)
            for f in database_files
        ]
        with ProgressBar():
            motifs = prune2df(databases, modules, str(motif_file), num_workers=workers)
        regulons = df2regulons(motifs)
        if not regulons:
            raise ValueError("No regulons retained after motif pruning.")
        motifs.to_csv(out_dir / "motifs.csv")
        with regulon_file.open("wb") as handle:
            pickle.dump(regulons, handle)

    # 4. Final targets and per-cell activity -----------------------------------
    if not regulons:
        raise ValueError("The regulon list is empty.")
    targets = [
        {"Regulon": r.name, "TF": r.transcription_factor,
         "Target": gene, "Weight": weight}
        for r in regulons for gene, weight in r.gene2weight.items()
    ]
    target_df = pd.DataFrame(targets)
    target_df.to_csv(out_dir / f"{cell_type}_final_regulon_targets.csv", index=False)
    target_df.groupby(["Regulon", "TF"])["Target"].nunique().rename(
        "n_targets").reset_index().to_csv(out_dir / "regulon_summary.csv", index=False)
    auc = aucell(ex_matrix, regulons, num_workers=workers)
    auc.to_csv(out_dir / f"{cell_type}_auc_matrix.csv")
    print(f"Saved {len(regulons)} regulons and AUC {auc.shape} to {out_dir}")


if __name__ == "__main__":
    main()
