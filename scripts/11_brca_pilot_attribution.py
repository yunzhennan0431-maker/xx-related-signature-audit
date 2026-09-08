"""
阶段2扩展:乳腺癌(BRCA)3个模型的基因-细胞类型归因。
数据复用前序课题"血管生成签名横向对比研究"泛癌种检验时已下载好的GSE176078
(Wu et al. 2021 BRCA scRNA-seq),不重复下载。方法与参考脚本
`/Users/yun/CRC_raw_data/pan_cancer/BRCA_GSE176078/phase2_brca_gse176078.py`一致
(log-normalize每个基因的counts,按celltype_major/minor分组求均值,跨细胞类型z-score)。
"""
import numpy as np
import pandas as pd
from scipy.io import mmread

DATA_DIR = "/Users/yun/CRC_raw_data/pan_cancer/BRCA_GSE176078/Wu_etal_2021_BRCA_scRNASeq/"
OUT_DIR = "analysis_output/data/"

MODELS = {
    "PMC13083393_autophagy_BRCA": ["ATG5", "SLC1A4", "AKT3", "CLTA", "HSPA2", "MAP2K5",
                                    "MAP1LC3B", "ITPR1", "DDIT4", "RPS6KA3", "DAPK1"],
    "PMC13406258_PCDgenes_BRCA": ["BRSK2", "CD24", "IFNG", "LAMB3", "PDX1", "PMAIP1", "SLC7A11", "TRIML2"],
    "PMC11855622_mitophagyE3_BRCA": ["ARIH1", "SIAH2", "UBR5", "WWP2"],
}
ALL_GENES = sorted(set(g for genes in MODELS.values() for g in genes))


def main():
    print("Loading genes/barcodes/metadata...")
    genes = pd.read_csv(DATA_DIR + "count_matrix_genes.tsv", header=None, names=["symbol"])["symbol"].tolist()
    barcodes = pd.read_csv(DATA_DIR + "count_matrix_barcodes.tsv", header=None, names=["barcode"])["barcode"].tolist()
    meta = pd.read_csv(DATA_DIR + "metadata.csv", index_col=0)
    print("meta shape:", meta.shape)
    print("celltype_major:", meta["celltype_major"].value_counts().to_dict())

    print("Loading sparse matrix (genes x cells)...")
    mat = mmread(DATA_DIR + "count_matrix_sparse.mtx").tocsr()
    print("matrix shape:", mat.shape)

    name_to_idx = {}
    for i, n in enumerate(genes):
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

    common = lognorm.index.intersection(meta.index)
    print("cells with matched metadata:", len(common), "/", len(barcodes))
    lognorm = lognorm.loc[common]
    lognorm["celltype_major"] = meta.loc[common, "celltype_major"].values

    mean_by_ct = lognorm.groupby("celltype_major")[present].mean().T
    n_by_ct = lognorm["celltype_major"].value_counts()
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
    out.to_csv(OUT_DIR + "phase2_brca_pilot_attribution.csv", index=False)
    print(f"\nn_cells by celltype_major:\n{n_by_ct}")
    print(f"\n写入 {OUT_DIR}phase2_brca_pilot_attribution.csv")


if __name__ == "__main__":
    main()
