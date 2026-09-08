"""
阶段2扩展:黑色素瘤(melanoma)、胶质瘤(glioma)、急性髓系白血病(AML)——
这3个癌种本地没有现成数据,从TISCH2(tisch.compbio.cn)新下载,存放于
/Users/yun/SignatureDL_raw_data/pan_cancer_phase2/(不进坚果云同步)。
方法与之前STAD/其他TISCH2格式数据集一致(`13_stad_pilot_attribution.py`同一套逻辑)。
"""
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from collections import Counter

OUT_DIR = "analysis_output/data/"
RAW_DIR = "/Users/yun/SignatureDL_raw_data/pan_cancer_phase2/"

DATASETS = {
    "SKCM_GSE120575": {
        "dir": RAW_DIR + "SKCM_GSE120575/",
        "models": {
            "PMC12219082_necroptosis_SKCM": ["TUFM", "CD53", "CLE2D", "KLRC1", "STAT4", "IFI35", "XCL1", "TAPBP", "SOD2"],
            "PMC12905907_disulfidptosis_SKCM": ["UCP2", "TMX4", "WWTR1", "CDH19", "CD79A"],
        },
    },
    "Glioma_GSE131928": {
        "dir": RAW_DIR + "Glioma_GSE131928/",
        "models": {
            "PMC9420977_cuproptosis_Glioma": ["FDX1", "LIPT1", "DLD", "DBT", "GCSH", "DLAT", "SLC31A1", "ATP7A", "ATP7B"],
            "PMC12239816_hypoxia_Glioma": ["ALOX5", "NGB", "SLC2A1", "ATF3", "OXSR1", "SLC2A3", "CAPG", "PRDX1",
                                            "SLC40A1", "CEBPG", "RGS4", "STEAP3", "HAMP", "RPL8", "TF", "HMOX1",
                                            "RRM2", "VEGFA", "MAP3K5", "SESN2", "XBP1", "NCF2", "SLC1A4"],
        },
    },
    "AML_GSE116256": {
        "dir": RAW_DIR + "AML_GSE116256/",
        "models": {
            "PMC12378281_disulfidptosis_AML": ["AIFM2", "FTMT", "MYB", "PRKAA2", "PSAT1", "SLC7A11"],
            "PMC12402452_exosome_AML": ["EXOSC4", "TMEM109", "THBS1", "MYH9", "HLA-DRA"],
        },
    },
}


def run_dataset(tag, cfg):
    data_dir = cfg["dir"]
    models = cfg["models"]
    print(f"\n########## {tag} ##########")
    print("Loading expression.h5...")
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
    return pd.concat(all_rows, ignore_index=True)


def main():
    all_dfs = []
    for tag, cfg in DATASETS.items():
        all_dfs.append(run_dataset(tag, cfg))
    out = pd.concat(all_dfs, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_newcancers_pilot_attribution.csv", index=False)
    print(f"\n写入 {OUT_DIR}phase2_newcancers_pilot_attribution.csv")


if __name__ == "__main__":
    main()
