"""
独立方法验证(响应"老板"反馈):现有验证锚点(21_validation_anchor_angiogenesis.py)只是
"回归一致性检验"——本研究z-score管线与前序课题z-score管线同源,不构成独立方法互证。

本脚本用一个数学上完全不同的统计算法重新做同一套归因,与z-score结果比对:
Wilcoxon秩和检验(one-vs-rest,per gene per celltype),取z统计量最大的细胞类型作为归因结果。
这正是Seurat FindAllMarkers、scanpy rank_genes_groups(method='wilcoxon')默认使用的marker基因
识别算法——单细胞领域最广泛使用、同行评议充分确立的标准方法,和z-score(细胞类型均值的
标准化)在数学原理上完全独立:前者是非参数秩检验,不假设分布形状,对离群值稳健;后者是
参数化的组间均值比较。

数据加载逻辑与21号脚本完全一致(同一批数据集、同一批benchmark基因),只是不在这一步把
表达矩阵坍缩成"按细胞类型求均值",而是保留逐细胞的表达值以计算秩。

基准集:前序血管生成课题13个模型(9个CRC模型45个不重复基因、4个泛癌种模型)已发表的
人工金标准,与21号脚本用的是同一个基准,确保两种方法的比较公平。
"""
import numpy as np
import pandas as pd
from scipy.io import mmread
import h5py
import scipy.sparse as sp

ANGIO_DIR = "/Users/yun/Library/CloudStorage/坚果云-yunzhennan0431@gmail.com/Nutstore/血管生成签名横向对比研究/analysis_output/data/"
OUT_DIR = "analysis_output/data/"
CRC_PAN = "/Users/yun/CRC_raw_data/pan_cancer/"
CRC81861 = "/Users/yun/CRC_raw_data/GSE81861/"
GSE178341_DIR = "/Users/yun/CRC_raw_data/GSE178341/"


def wilcoxon_attribution(expr_df: pd.DataFrame, celltype: pd.Series, genes) -> pd.DataFrame:
    """对expr_df中的每个基因,逐细胞类型做one-vs-rest Wilcoxon秩和检验(含结列修正的
    正态近似z统计量),取z最大的细胞类型作为归因结果。数学上等价于scanpy/Seurat默认的
    marker基因打分方法,与z-score-of-means是完全不同的统计路线。"""
    n_total = len(celltype)
    cts = celltype.unique()
    ct_values = celltype.values
    rows = []
    for g in genes:
        if g not in expr_df.columns:
            continue
        vals = expr_df[g].values.astype(float)
        ranks = pd.Series(vals).rank(method="average").values
        _, tie_counts = np.unique(ranks, return_counts=True)
        tie_term = np.sum(tie_counts ** 3 - tie_counts)
        best_ct, best_z = None, -np.inf
        for ct in cts:
            mask = ct_values == ct
            n1 = int(mask.sum())
            n2 = n_total - n1
            if n1 == 0 or n2 == 0:
                continue
            R1 = ranks[mask].sum()
            U1 = R1 - n1 * (n1 + 1) / 2.0
            mu = n1 * n2 / 2.0
            if n_total > 1:
                sigma_sq = (n1 * n2 / 12.0) * ((n_total + 1) - tie_term / (n_total * (n_total - 1)))
            else:
                sigma_sq = 0
            sigma = np.sqrt(sigma_sq) if sigma_sq > 0 else 0
            z = (U1 - mu) / sigma if sigma > 0 else 0.0
            if z > best_z:
                best_z, best_ct = z, ct
        rows.append({"gene": g, "top_celltype": best_ct, "wilcoxon_z": round(best_z, 2)})
    return pd.DataFrame(rows).set_index("gene")


def load_gse81861_cells(genes):
    def load(fp):
        df = pd.read_csv(fp, index_col=0)
        df.index = [str(i).split("_")[1] if len(str(i).split("_")) >= 3 else str(i) for i in df.index]
        return df

    nm = load(CRC81861 + "GSE81861_CRC_NM_all_cells_FPKM.csv")
    tm = load(CRC81861 + "GSE81861_CRC_tumor_all_cells_FPKM.csv")
    nm = nm.groupby(nm.index).mean()
    tm = tm.groupby(tm.index).mean()
    common_genes = nm.index.intersection(tm.index)
    nm, tm = nm.loc[common_genes], tm.loc[common_genes]
    combined = pd.concat([nm, tm], axis=1)
    cell_type = [c.split("__")[1] if "__" in c else "NA" for c in combined.columns]
    meta = pd.DataFrame({"cell": combined.columns, "celltype": cell_type})
    meta = meta[meta.celltype != "NA"]
    combined = combined[meta.cell.values]
    genes_present = [g for g in genes if g in combined.index]
    logexpr = np.log2(combined.loc[genes_present] + 1).T
    logexpr.index = meta.cell.values
    celltype = pd.Series(meta.celltype.values, index=meta.cell.values)
    return logexpr, celltype


def load_gse178341_cells(genes):
    f = h5py.File(GSE178341_DIR + "GSE178341_crc10x_full_c295v4_submit.h5", "r")
    n_genes, n_cells = f["matrix/shape"][:]
    names = np.array([x.decode() for x in f["matrix/features/name"][:]])
    barcodes = np.array([x.decode() for x in f["matrix/barcodes"][:]])
    name_to_idx = {}
    for i, n in enumerate(names):
        name_to_idx.setdefault(n, []).append(i)
    present = [g for g in genes if g in name_to_idx]

    data = f["matrix/data"][:]; indices = f["matrix/indices"][:]; indptr = f["matrix/indptr"][:]
    mat = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).tocsr()
    total_counts = np.asarray(mat.sum(axis=0)).flatten()
    total_counts[total_counts == 0] = 1
    gene_expr = {}
    for g in present:
        idxs = name_to_idx[g]
        gene_expr[g] = np.asarray(mat[idxs, :].sum(axis=0)).flatten()
    expr_df = pd.DataFrame(gene_expr, index=barcodes)
    norm = expr_df.div(total_counts, axis=0) * 1e4
    lognorm = np.log1p(norm)

    cluster = pd.read_csv(GSE178341_DIR + "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz", index_col=0)
    common = lognorm.index.intersection(cluster.index)
    lognorm = lognorm.loc[common]
    celltype = pd.Series(cluster.loc[common, "clMidwayPr"].values, index=common)
    return lognorm, celltype


def load_brca_cells(genes):
    data_dir = CRC_PAN + "BRCA_GSE176078/Wu_etal_2021_BRCA_scRNASeq/"
    gene_list = pd.read_csv(data_dir + "count_matrix_genes.tsv", header=None, names=["symbol"])["symbol"].tolist()
    barcodes = pd.read_csv(data_dir + "count_matrix_barcodes.tsv", header=None, names=["barcode"])["barcode"].tolist()
    meta = pd.read_csv(data_dir + "metadata.csv", index_col=0)
    mat = mmread(data_dir + "count_matrix_sparse.mtx").tocsr()
    name_to_idx = {}
    for i, n in enumerate(gene_list):
        name_to_idx.setdefault(n, []).append(i)
    present = [g for g in genes if g in name_to_idx]
    total_counts = np.asarray(mat.sum(axis=0)).flatten()
    total_counts[total_counts == 0] = 1
    gene_expr = {}
    for g in present:
        idxs = name_to_idx[g]
        gene_expr[g] = np.asarray(mat[idxs, :].sum(axis=0)).flatten()
    expr_df = pd.DataFrame(gene_expr, index=barcodes)
    norm = expr_df.div(total_counts, axis=0) * 1e4
    lognorm = np.log1p(norm)
    common = lognorm.index.intersection(meta.index)
    lognorm = lognorm.loc[common]
    celltype = pd.Series(meta.loc[common, "celltype_major"].values, index=common)
    return lognorm, celltype


def load_hcc_cells(genes):
    data_dir = CRC_PAN + "HCC_GSE151530/"
    gene_tab = pd.read_csv(data_dir + "GSE151530_genes.tsv.gz", sep="\t", header=None, names=["ensembl", "symbol"])
    barcodes = pd.read_csv(data_dir + "GSE151530_barcodes.tsv.gz", header=None, names=["barcode"])["barcode"].tolist()
    info = pd.read_csv(data_dir + "info.txt", sep="\t").set_index("Cell")
    mat = mmread(data_dir + "GSE151530_matrix.mtx.gz").tocsr()
    name_to_idx = {}
    for i, n in enumerate(gene_tab["symbol"]):
        name_to_idx.setdefault(n, []).append(i)
    present = [g for g in genes if g in name_to_idx]
    gene_expr = {}
    for g in present:
        idxs = name_to_idx[g]
        gene_expr[g] = np.asarray(mat[idxs, :].sum(axis=0)).flatten()
    expr_df = pd.DataFrame(gene_expr, index=barcodes)
    total_counts = np.asarray(mat.sum(axis=0)).flatten()
    total_counts[total_counts == 0] = 1
    norm = expr_df.div(total_counts, axis=0) * 1e4
    lognorm = np.log1p(norm)
    common = lognorm.index.intersection(info.index)
    lognorm = lognorm.loc[common]
    celltype = pd.Series(info.loc[common, "Type"].values, index=common)
    keep = celltype != "unclassified"
    return lognorm.loc[keep], celltype.loc[keep]


def load_tisch2_cells(data_dir, genes):
    f = h5py.File(data_dir + "expression.h5", "r")
    n_genes, n_cells = f["matrix/shape"][:]
    names = np.array([x.decode() for x in f["matrix/features/name"][:]])
    barcodes = np.array([x.decode() for x in f["matrix/barcodes"][:]])
    name_to_idx = {}
    for i, n in enumerate(names):
        name_to_idx.setdefault(n, []).append(i)
    present = [g for g in genes if g in name_to_idx]
    data = f["matrix/data"][:]; indices = f["matrix/indices"][:]; indptr = f["matrix/indptr"][:]
    mat = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).tocsr()
    gene_expr = {}
    for g in present:
        idxs = name_to_idx[g]
        gene_expr[g] = np.asarray(mat[idxs, :].max(axis=0).todense()).flatten()
    expr_df = pd.DataFrame(gene_expr, index=barcodes)
    meta = pd.read_csv(data_dir + "meta.tsv", sep="\t").set_index("Cell")
    common = expr_df.index.intersection(meta.index)
    expr_df = expr_df.loc[common]
    celltype = pd.Series(meta.loc[common, "Celltype (major-lineage)"].values, index=common)
    return expr_df, celltype


def compare(new_df, gold_path, zscore_path, label):
    gold = pd.read_csv(gold_path, index_col=0)
    common = new_df.index.intersection(gold.index)
    merged = pd.DataFrame({
        "gene": common,
        "wilcoxon_top_celltype": new_df.loc[common, "top_celltype"].values,
        "gold_top_celltype": gold.loc[common, "top_celltype"].values,
    })
    merged["match_vs_gold"] = merged["wilcoxon_top_celltype"] == merged["gold_top_celltype"]
    n_match = merged["match_vs_gold"].sum()
    print(f"\n=== {label}: Wilcoxon vs 前序课题金标准 {n_match}/{len(merged)} 一致 ({n_match/len(merged)*100:.1f}%) ===")
    print(merged.to_string(index=False))
    return merged


def main():
    model_df = pd.read_csv(ANGIO_DIR + "model_gene_lists.csv")
    crc_genes = sorted(model_df["gene"].unique())
    pancancer_df = pd.read_csv(ANGIO_DIR + "pancancer_model_gene_lists.csv")

    all_results = []

    print("########## GSE81861 (CRC, 45基因, Wilcoxon独立归因) ##########")
    expr, ct = load_gse81861_cells(crc_genes)
    new_81861 = wilcoxon_attribution(expr, ct, crc_genes)
    m = compare(new_81861, ANGIO_DIR + "phase1_gse81861_gene_celltype_attribution.csv",
                None, "GSE81861")
    m["dataset"] = "GSE81861"; all_results.append(m)

    print("\n########## GSE178341 clMidwayPr (CRC, 45基因, Wilcoxon独立归因) ##########")
    expr, ct = load_gse178341_cells(crc_genes)
    new_178341 = wilcoxon_attribution(expr, ct, crc_genes)
    m = compare(new_178341, ANGIO_DIR + "phase1_gse178341_gene_midway_attribution.csv",
                None, "GSE178341")
    m["dataset"] = "GSE178341"; all_results.append(m)

    print("\n########## KDR跨数据集一致性检查(前序课题核心发现,独立算法复现) ##########")
    kdr_81861 = new_81861.loc["KDR", "top_celltype"] if "KDR" in new_81861.index else None
    kdr_178341 = new_178341.loc["KDR", "top_celltype"] if "KDR" in new_178341.index else None
    print(f"KDR in GSE81861 (Wilcoxon): {kdr_81861}; KDR in GSE178341 clMidwayPr (Wilcoxon): {kdr_178341}")
    endo_labels = {"Endothelial", "Endo"}
    print(f"KDR两个数据集均落在内皮细胞(Wilcoxon独立算法): {kdr_81861 in endo_labels and kdr_178341 in endo_labels}")

    for cancer, genes, loader, gold_file in [
        ("BRCA", pancancer_df[pancancer_df.cancer_type.str.contains("BRCA")]["gene"].tolist(),
         load_brca_cells, "pancancer_brca_gse176078_major_attribution.csv"),
        ("HCC", pancancer_df[pancancer_df.cancer_type == "HCC"]["gene"].tolist(),
         load_hcc_cells, "pancancer_hcc_gse151530_attribution.csv"),
        ("STAD", pancancer_df[pancancer_df.cancer_type.str.contains("STAD")]["gene"].tolist(),
         lambda g: load_tisch2_cells(CRC_PAN + "STAD_GSE167297/", g), "pancancer_stad_gse167297_attribution.csv"),
        ("BLCA", pancancer_df[pancancer_df.cancer_type.str.contains("BLCA")]["gene"].tolist(),
         lambda g: load_tisch2_cells(CRC_PAN + "BLCA_GSE130001/", g), "pancancer_blca_gse130001_attribution.csv"),
    ]:
        print(f"\n########## {cancer}(泛癌种模型,{len(genes)}基因, Wilcoxon独立归因) ##########")
        expr, ct = loader(genes)
        new_df = wilcoxon_attribution(expr, ct, genes)
        m = compare(new_df, ANGIO_DIR + gold_file, None, cancer)
        m["dataset"] = cancer; all_results.append(m)

    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(OUT_DIR + "independent_wilcoxon_validation_comparison.csv", index=False)
    overall = combined.dropna(subset=["wilcoxon_top_celltype", "gold_top_celltype"])
    n_match = overall["match_vs_gold"].sum()
    print(f"\n\n========== Wilcoxon独立算法 vs 前序课题金标准 总体一致率: {n_match}/{len(overall)} ({n_match/len(overall)*100:.1f}%) ==========")
    print("\n按数据集分组一致率:")
    print(combined.groupby("dataset")["match_vs_gold"].agg(["sum", "count"]))

    # 与本研究z-score管线结果的一致率(方法互证的核心比较)
    zscore_cmp = pd.read_csv(OUT_DIR + "validation_anchor_angiogenesis_comparison.csv")
    zscore_cmp = zscore_cmp.rename(columns={"new_top_celltype": "zscore_top_celltype"})
    merge2 = combined.merge(zscore_cmp[["gene", "dataset", "zscore_top_celltype"]], on=["gene", "dataset"], how="inner")
    merge2["wilcoxon_vs_zscore_match"] = merge2["wilcoxon_top_celltype"] == merge2["zscore_top_celltype"]
    n_wz = merge2["wilcoxon_vs_zscore_match"].sum()
    print(f"\n========== Wilcoxon独立算法 vs 本研究z-score管线 总体一致率: {n_wz}/{len(merge2)} ({n_wz/len(merge2)*100:.1f}%) ==========")
    print(merge2.groupby("dataset")["wilcoxon_vs_zscore_match"].agg(["sum", "count"]))
    merge2.to_csv(OUT_DIR + "independent_wilcoxon_vs_zscore_comparison.csv", index=False)

    print(f"\n写入 {OUT_DIR}independent_wilcoxon_validation_comparison.csv")
    print(f"写入 {OUT_DIR}independent_wilcoxon_vs_zscore_comparison.csv")


if __name__ == "__main__":
    main()
