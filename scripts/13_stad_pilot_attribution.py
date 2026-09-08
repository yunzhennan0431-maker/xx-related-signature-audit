"""
阶段2扩展:胃癌(STAD)模型的基因-细胞类型归因。
数据复用前序课题泛癌种检验时已下载好的GSE167297(TISCH2格式),方法与
`/Users/yun/CRC_raw_data/pan_cancer/tisch2_attribution_template.py`通用模板一致。
lncRNA基因名注意点号/短横线两种写法都试(XML导出常把'-'转成'.').
"""
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

DATA_DIR = "/Users/yun/CRC_raw_data/pan_cancer/STAD_GSE167297/"
OUT_DIR = "analysis_output/data/"

MODELS = {
    "PMC11760997_senescence_STAD": ["LINC01579", "LINC02326", "C3orf36", "CACNA1C.AS4", "DNPEP.AS1", "FCHO2.DT",
                                     "GREM1.AS1", "HHIP.AS1", "LINC01344", "LINC01589", "LINC01615", "LINC02864",
                                     "LINC02975", "MAGEA6.DT", "STARD4.AS1", "TICAM2.AS1"],
    "PMC12399450_mitochondrial_STAD": ["ATP8A2", "TARS2", "COX15"],
}


def main():
    print("Loading expression.h5...")
    f = h5py.File(DATA_DIR + "expression.h5", "r")
    n_genes, n_cells = f["matrix/shape"][:]
    names = np.array([x.decode() for x in f["matrix/features/name"][:]])
    barcodes = np.array([x.decode() for x in f["matrix/barcodes"][:]])
    print("genes:", n_genes, "cells:", n_cells)

    name_to_idx = {}
    for i, n in enumerate(names):
        name_to_idx.setdefault(n, []).append(i)

    # 点号/短横线两种写法都试一遍,取能匹配上的
    all_raw_genes = sorted(set(g for genes in MODELS.values() for g in genes))
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
    meta = pd.read_csv(DATA_DIR + "meta.tsv", sep="\t")
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

    from collections import Counter
    all_rows = []
    for model_name, raw_genes in MODELS.items():
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

    out = pd.concat(all_rows, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_stad_pilot_attribution.csv", index=False)
    print(f"\nn_cells by type:\n{n_by_ct}")
    print(f"\n写入 {OUT_DIR}phase2_stad_pilot_attribution.csv")


if __name__ == "__main__":
    main()
