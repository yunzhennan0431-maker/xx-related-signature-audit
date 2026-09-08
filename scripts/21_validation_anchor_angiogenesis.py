"""
阶段2归因方法第一层验证锚点:用本课题当前的z-score归因管线(独立实现,和批次6-8同一套
方法),重新跑一遍前序课题("血管生成签名横向对比研究")已发表、经人工核实过的13个血管
生成相关模型(9个CRC模型+4个泛癌种模型),把逐基因归因结果与前序课题已发表的金标准结果
逐一比对,报告一致率。这是`跨主题meta分析_设计方案.md`/`论文初稿_中文.md`第2.4节里明确
标注"尚未执行"的工作,本脚本补上这一步。

前序课题原始基因列表/金标准归因结果只读取,不修改、不重新分发。
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
SDL = "/Users/yun/SignatureDL_raw_data/pan_cancer_phase2/"


def zscore_attribution(mean_by_ct: pd.DataFrame) -> pd.DataFrame:
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)
    return pd.DataFrame({"top_celltype": top_ct, "zscore": z.max(axis=1).round(2)})


def run_gse81861(genes):
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
    logexpr = np.log2(combined.loc[genes_present] + 1)
    result = {}
    for ct in meta.celltype.unique():
        cells = meta[meta.celltype == ct].cell.values
        result[ct] = logexpr[cells].mean(axis=1)
    mean_by_ct = pd.DataFrame(result)
    out = zscore_attribution(mean_by_ct)
    return out.reindex(genes)


def run_gse178341(genes):
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
    lognorm["clMidwayPr"] = cluster.loc[common, "clMidwayPr"].values

    mean_by_ct = lognorm.groupby("clMidwayPr")[present].mean().T
    out = zscore_attribution(mean_by_ct)
    return out.reindex(genes)


def run_brca(genes):
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
    lognorm["celltype"] = meta.loc[common, "celltype_major"].values
    mean_by_ct = lognorm.groupby("celltype")[present].mean().T
    out = zscore_attribution(mean_by_ct)
    return out.reindex(genes)


def run_hcc(genes):
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
    lognorm["celltype"] = info.loc[common, "Type"].values
    lognorm = lognorm[lognorm["celltype"] != "unclassified"]
    mean_by_ct = lognorm.groupby("celltype")[present].mean().T
    out = zscore_attribution(mean_by_ct)
    return out.reindex(genes)


def run_tisch2(data_dir, genes):
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
    expr_df["celltype"] = meta.loc[common, "Celltype (major-lineage)"].values
    mean_by_ct = expr_df.groupby("celltype")[present].mean().T
    out = zscore_attribution(mean_by_ct)
    return out.reindex(genes)


def compare(new_df, gold_path, label):
    gold = pd.read_csv(gold_path, index_col=0)
    common = new_df.index.intersection(gold.index)
    merged = pd.DataFrame({
        "gene": common,
        "new_top_celltype": new_df.loc[common, "top_celltype"].values,
        "gold_top_celltype": gold.loc[common, "top_celltype"].values,
    })
    merged["match"] = merged["new_top_celltype"] == merged["gold_top_celltype"]
    n_match = merged["match"].sum()
    print(f"\n=== {label}: {n_match}/{len(merged)} 一致 ({n_match/len(merged)*100:.1f}%) ===")
    print(merged.to_string(index=False))
    return merged


def main():
    model_df = pd.read_csv(ANGIO_DIR + "model_gene_lists.csv")
    crc_genes = sorted(model_df["gene"].unique())
    pancancer_df = pd.read_csv(ANGIO_DIR + "pancancer_model_gene_lists.csv")

    all_results = []

    print("########## GSE81861 (CRC, 9模型/45基因) ##########")
    new_81861 = run_gse81861(crc_genes)
    m = compare(new_81861, ANGIO_DIR + "phase1_gse81861_gene_celltype_attribution.csv", "GSE81861")
    m["dataset"] = "GSE81861"; all_results.append(m)

    print("\n########## GSE178341 clMidwayPr (CRC, 45基因) ##########")
    new_178341 = run_gse178341(crc_genes)
    m = compare(new_178341, ANGIO_DIR + "phase1_gse178341_gene_midway_attribution.csv", "GSE178341")
    m["dataset"] = "GSE178341"; all_results.append(m)

    print("\n########## KDR跨数据集一致性检查(前序课题核心发现的复现) ##########")
    kdr_81861 = new_81861.loc["KDR", "top_celltype"] if "KDR" in new_81861.index else None
    kdr_178341 = new_178341.loc["KDR", "top_celltype"] if "KDR" in new_178341.index else None
    print(f"KDR in GSE81861: {kdr_81861}; KDR in GSE178341(clMidwayPr): {kdr_178341}")
    endo_labels = {"Endothelial", "Endo"}
    print(f"KDR两个数据集均落在内皮细胞: {kdr_81861 in endo_labels and kdr_178341 in endo_labels}")

    for cancer, genes, runner, gold_file in [
        ("BRCA", pancancer_df[pancancer_df.cancer_type.str.contains("BRCA")]["gene"].tolist(),
         lambda g: run_brca(g), "pancancer_brca_gse176078_major_attribution.csv"),
        ("HCC", pancancer_df[pancancer_df.cancer_type == "HCC"]["gene"].tolist(),
         lambda g: run_hcc(g), "pancancer_hcc_gse151530_attribution.csv"),
        ("STAD", pancancer_df[pancancer_df.cancer_type.str.contains("STAD")]["gene"].tolist(),
         lambda g: run_tisch2(CRC_PAN + "STAD_GSE167297/", g), "pancancer_stad_gse167297_attribution.csv"),
        ("BLCA", pancancer_df[pancancer_df.cancer_type.str.contains("BLCA")]["gene"].tolist(),
         lambda g: run_tisch2(CRC_PAN + "BLCA_GSE130001/", g), "pancancer_blca_gse130001_attribution.csv"),
    ]:
        print(f"\n########## {cancer}(泛癌种模型,{len(genes)}基因) ##########")
        new_df = runner(genes)
        m = compare(new_df, ANGIO_DIR + gold_file, cancer)
        m["dataset"] = cancer; all_results.append(m)

    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(OUT_DIR + "validation_anchor_angiogenesis_comparison.csv", index=False)
    overall = combined.dropna(subset=["new_top_celltype", "gold_top_celltype"])
    n_match = overall["match"].sum()
    print(f"\n\n========== 总体一致率: {n_match}/{len(overall)} ({n_match/len(overall)*100:.1f}%) ==========")
    print("\n按数据集分组一致率:")
    print(combined.groupby("dataset")["match"].agg(["sum", "count"]))
    print(f"\n写入 {OUT_DIR}validation_anchor_angiogenesis_comparison.csv")


if __name__ == "__main__":
    main()
