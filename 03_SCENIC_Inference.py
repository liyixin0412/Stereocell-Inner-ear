"""Representative pySCENIC workflow for SGN, Glia or Epi.
Input: h5ad with raw counts in layers['counts']; X is the input used for
the original filter/normalize/log1p/HVG procedure and must be unlogged counts.
Inference parameters are retained; worker count is configurable.
"""
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-type", choices=["SGN", "Glia", "Epi"], default="SGN")
    parser.add_argument("--input")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    import os, glob, pickle
    import numpy as np
    import pandas as pd
    import scipy.sparse
    import scanpy as sc
    from dask.distributed import Client
    from dask.diagnostics import ProgressBar
    from arboreto.utils import load_tf_names
    from arboreto.algo import grnboost2
    from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
    from pyscenic.utils import modules_from_adjacencies
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell

    cell_type = args.cell_type
    adata_path = args.input or f"data/Merge-{cell_type}-Final-pySCENIC.h5ad"
    DATABASES_GLOB = "references/mm10__refseq-r80__*tss.mc9nr.genes_vs_motifs.rankings.feather"
    MOTIF_ANNOTATIONS_FNAME = "references/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
    tf_list_path = "references/mm_mgi_tfs.txt"
    OUTPUT_DIR = f"results/03/{cell_type}"
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    MOTIFS_FNAME = os.path.join(OUTPUT_DIR, "motifs.csv")
    REGULONS_FNAME = os.path.join(OUTPUT_DIR, "regulons.p")
    adata = sc.read_h5ad(adata_path)
    counts = adata.layers["counts"]
    ex_matrix = pd.DataFrame(counts.toarray() if scipy.sparse.issparse(counts) else counts,
                             index=adata.obs_names, columns=adata.var_names)
    tf_names = load_tf_names(tf_list_path)


    import scanpy as sc



    sc.pp.filter_genes(adata, min_cells=10)  # 过滤在少于10个细胞中表达的基因
    sc.pp.normalize_total(adata, target_sum=1e4)  # 标准化
    sc.pp.log1p(adata)  # log转换

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=5000,  # 选择5000个高变基因
        flavor='seurat',   # 使用Seurat的方法
        batch_key=None     # 如果有批次效应，可以指定批次列
    )


    np.random.seed(42)
    hvg_genes = adata.var_names[adata.var.highly_variable]
    ex_matrix_hvg = ex_matrix.loc[:, ex_matrix.columns.intersection(hvg_genes)]


    client = Client(n_workers=args.workers, threads_per_worker=1)
    adjacencies = grnboost2(
        expression_data=ex_matrix_hvg,
        tf_names=tf_names,
        client_or_address=client,
        verbose=True# 显式传递 client
    )
    adjacencies.to_csv(os.path.join(OUTPUT_DIR, f"{cell_type}_adjacencies.csv"), index=False)
    client.close()


    db_fnames = sorted(glob.glob(DATABASES_GLOB))
    if not db_fnames:
        raise FileNotFoundError(DATABASES_GLOB)
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
    dbs

    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))


    with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, num_workers=args.workers)

    regulons = df2regulons(df)
    df.to_csv(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "wb") as f:
        pickle.dump(regulons, f)


    auc_mtx = aucell(ex_matrix, regulons, num_workers=args.workers)
    auc_mtx.to_csv(os.path.join(OUTPUT_DIR, f"{cell_type}_auc_matrix.csv"))

    # Export final motif-pruned targets, not raw adjacency edges.
    targets = [{"Regulon": r.name, "TF": r.transcription_factor,
                "Target": gene, "Weight": weight}
               for r in regulons for gene, weight in r.gene2weight.items()]
    pd.DataFrame(targets).to_csv(
        os.path.join(OUTPUT_DIR, f"{cell_type}_regulon_targets.csv"), index=False)


if __name__ == "__main__":
    main()
