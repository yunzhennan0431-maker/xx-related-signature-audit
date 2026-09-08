"""
阶段2继续扩大:肾癌(KIRC)、骨肉瘤(OS)、子宫内膜癌(UCEC)、胰腺癌(PAAD)、
肺腺癌(NSCLC)、宫颈癌(CESC)——6个此前只有1个模型、本地没有数据的癌种,
一次性从TISCH2补齐,每个癌种凑够至少1个独立数据点。
数据存于/Users/yun/SignatureDL_raw_data/pan_cancer_phase2/,不进坚果云。
"""
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from collections import Counter

OUT_DIR = "analysis_output/data/"
RAW_DIR = "/Users/yun/SignatureDL_raw_data/pan_cancer_phase2/"

DATASETS = {
    "KIRC_GSE159115": {
        "dir": RAW_DIR + "KIRC_GSE159115/",
        "models": {"PMC10184413_cuproptosis_ferroptosis_KIRC": ["TRIB3", "SLC2A3", "PML", "CD44", "CDKN2A", "MIOX"]},
    },
    "OS_GSE162454": {
        "dir": RAW_DIR + "OS_GSE162454/",
        "models": {"PMC12047218_glycolysis_OS": ["CHPF", "RRAGD", "TPR", "VCAN"]},
    },
    "UCEC_GSE139555": {
        "dir": RAW_DIR + "UCEC_GSE139555/",
        "models": {"PMC10469729_glutamine_UCEC": ["PHGDH", "OTC", "ASRGL1", "ASNS", "NR1H4"]},
    },
    "PAAD_GSE154778": {
        "dir": RAW_DIR + "PAAD_GSE154778/",
        "models": {"PMC11504279_hypoxia_PAAD": ["SLC2A1", "CA9", "PLAU"]},
    },
    "NSCLC_GSE131907": {
        "dir": RAW_DIR + "NSCLC_GSE131907/",
        "models": {"PMC13333213_pyroptosis_NSCLC": ["PRKACA", "CASP4"]},
    },
    "CESC_GSE168652": {
        "dir": RAW_DIR + "CESC_GSE168652/",
        "models": {"PMC10564413_senescence_CESC": ["PTTG1", "E2F1", "CBX7", "SPOP", "RPS6KA6", "ABI3"]},
    },
}


def run_dataset(tag, cfg):
    data_dir = cfg["dir"]
    models = cfg["models"]
    print(f"\n########## {tag} ##########")
    f = h5py.File(data_dir + "expression.h5", "r")
    n_genes, n_cells = f["matrix/shape"][:]
    names = np.array([x.decode() for x in f["matrix/features/name"][:]])
    barcodes = np.array([x.decode() for x in f["matrix/barcodes"][:]])
    print("genes:", n_genes, "cells:", n_cells)

    name_to_idx = {}
    for i, n in enumerate(names):
        name_to_idx.setdefault(n, []).append(i)

    all_raw_genes = sorted(set(g for genes in models.values() for g in genes))
    resolved = {}
    for g in all_raw_genes:
        alt = g.replace(".", "-")
        if g in name_to_idx:
            resolved[g] = g
        elif alt in name_to_idx:
            resolved[g] = alt
    present = list(resolved.values())
    print("present:", len(present), "/", len(all_raw_genes))

    data = f["matrix/data"][:]
    indices = f["matrix/indices"][:]
    indptr = f["matrix/indptr"][:]
    mat = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells))
    del data, indices, indptr
    mat_csr = mat.tocsr()

    gene_expr = {}
    for g in present:
        idxs = name_to_idx[g]
        row = np.asarray(mat_csr[idxs, :].max(axis=0).todense()).flatten()
        gene_expr[g] = row
    del mat, mat_csr

    expr_df = pd.DataFrame(gene_expr, index=barcodes)
    meta = pd.read_csv(data_dir + "meta.tsv", sep="\t")
    meta = meta.set_index("Cell")

    common = expr_df.index.intersection(meta.index)
    print("cells with matched metadata:", len(common), "/", len(barcodes))
    expr_df = expr_df.loc[common]
    celltype = meta.loc[common, "Celltype (major-lineage)"]
    expr_df["celltype"] = celltype.values

    mean_by_ct = expr_df.groupby("celltype")[present].mean().T
    n_by_ct = expr_df["celltype"].value_counts()
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)

    all_rows = []
    for model_name, raw_genes in models.items():
        model_present = [resolved[g] for g in raw_genes if g in resolved]
        print(f"\n=== {model_name} ({len(model_present)}/{len(raw_genes)}基因在数据集中找到) ===")
        sub = pd.DataFrame({
            "model": model_name, "gene": model_present,
            "top_celltype": top_ct[model_present].values,
            "zscore": z.loc[model_present].max(axis=1).round(2).values,
        })
        all_rows.append(sub)
        print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
        print("细胞类型归因分布:", dict(Counter(sub.top_celltype)))

    print(f"\nn_cells by type:\n{n_by_ct}")
    print(f"\n[DATASET_CELLTYPE_COUNTS entry for {tag}]:")
    print(dict(n_by_ct))
    return pd.concat(all_rows, ignore_index=True)


def main():
    all_dfs = []
    for tag, cfg in DATASETS.items():
        all_dfs.append(run_dataset(tag, cfg))
    out = pd.concat(all_dfs, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_morecancers_pilot_attribution.csv", index=False)
    print(f"\n写入 {OUT_DIR}phase2_morecancers_pilot_attribution.csv")


if __name__ == "__main__":
    main()
