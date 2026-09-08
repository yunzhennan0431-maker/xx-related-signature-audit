"""
阶段2扩展:肝癌(HCC)2个模型的基因-细胞类型归因。
数据复用前序课题泛癌种检验时已下载好的GSE151530,方法与参考脚本
`/Users/yun/CRC_raw_data/pan_cancer/HCC_GSE151530/phase2_hcc_gse151530.py`一致。
"""
import numpy as np
import pandas as pd
from scipy.io import mmread

DATA_DIR = "/Users/yun/CRC_raw_data/pan_cancer/HCC_GSE151530/"
OUT_DIR = "analysis_output/data/"

MODELS = {
    "PMC12647489_cholesterol_HCC": ["EPHX2", "FABP5", "SQLE", "ANXA2", "ADH4", "HMGCS2", "CYP7A1", "ACADL"],
    "PMC13490468_ATIC_HCC": ["ATIC", "RHEB", "TMEM74", "PRKCD"],
    "PMC9906058_ferroptosis_cuproptosis_HCC": ["TXNRD1", "FTL", "GPX4", "PRDX1", "VDAC2", "OTUB1", "NRAS", "SLC38A1", "SLC1A5"],
    "PMC11577587_glutamine_HCC": ["G6PD", "GPX7", "RRM1", "RRM2"],
}
ALL_GENES = sorted(set(g for genes in MODELS.values() for g in genes))


def main():
    print("Loading genes/barcodes...")
    genes = pd.read_csv(DATA_DIR + "GSE151530_genes.tsv.gz", sep="\t", header=None, names=["ensembl", "symbol"])
    barcodes = pd.read_csv(DATA_DIR + "GSE151530_barcodes.tsv.gz", header=None, names=["barcode"])["barcode"].tolist()

    info = pd.read_csv(DATA_DIR + "info.txt", sep="\t")
    info = info.set_index("Cell")
    print("info shape:", info.shape, "cell types:", info["Type"].value_counts().to_dict())

    print("Loading sparse matrix (genes x cells)...")
    mat = mmread(DATA_DIR + "GSE151530_matrix.mtx.gz").tocsr()
    print("matrix shape:", mat.shape)

    name_to_idx = {}
    for i, n in enumerate(genes["symbol"]):
        name_to_idx.setdefault(n, []).append(i)

    present = [g for g in ALL_GENES if g in name_to_idx]
    missing = [g for g in ALL_GENES if g not in name_to_idx]
    print("present:", len(present), "missing:", missing)

    total_counts = np.asarray(mat.sum(axis=0)).flatten()
    total_counts[total_counts == 0] = 1

    gene_expr = {}
    for g in present:
        idxs = name_to_idx[g]
        row = np.asarray(mat[idxs, :].sum(axis=0)).flatten()
        gene_expr[g] = row

    expr_df = pd.DataFrame(gene_expr, index=barcodes)
    norm = expr_df.div(total_counts, axis=0) * 1e4
    lognorm = np.log1p(norm)

    common = lognorm.index.intersection(info.index)
    print("cells with matched cell type:", len(common), "/", len(barcodes))
    lognorm = lognorm.loc[common]
    celltype = info.loc[common, "Type"]
    lognorm["celltype"] = celltype.values
    lognorm = lognorm[lognorm["celltype"] != "unclassified"]

    mean_by_ct = lognorm.groupby("celltype")[present].mean().T
    n_by_ct = lognorm["celltype"].value_counts()
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)

    all_rows = []
    for model_name, model_genes in MODELS.items():
        model_genes_present = [g for g in model_genes if g in present]
        print(f"\n=== {model_name} ({len(model_genes_present)}/{len(model_genes)}基因在数据集中找到) ===")
        sub = pd.DataFrame({
            "model": model_name, "gene": model_genes_present,
            "top_celltype": top_ct[model_genes_present].values,
            "zscore": z.loc[model_genes_present].max(axis=1).round(2).values,
        })
        all_rows.append(sub)
        print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
        from collections import Counter
        print("细胞类型归因分布:", dict(Counter(sub.top_celltype)))

    out = pd.concat(all_rows, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_hcc_pilot_attribution.csv", index=False)
    print(f"\nn_cells by type (excl. unclassified):\n{n_by_ct}")
    print(f"\n写入 {OUT_DIR}phase2_hcc_pilot_attribution.csv")


if __name__ == "__main__":
    main()
