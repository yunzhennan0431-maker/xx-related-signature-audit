"""
用独立算法(Wilcoxon秩和检验,见40号脚本)重新对全部108个模型、全部21个数据集做归因,
检验论文headline结论(23.1%显著率、H1a:H1b比值)对"归因算法选择"是否稳健——这是40号脚本
只在119个验证锚点基准基因上做独立算法比对之后,更直接回答"会不会动摇论文结论本身"的问题。

复用15号脚本(跨主题meta分析)的MODEL_TO_DATASET/DATASET_CELLTYPE_COUNTS/CELLTYPE_CATEGORY/
TAUTOLOGICAL_MODELS和analyze()函数不做任何修改,只把输入的归因结果来源从z-score换成Wilcoxon,
确保统计口径完全一致、比较公平。
"""
import importlib.util
import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sp
from scipy.io import mmread

spec15 = importlib.util.spec_from_file_location("m15", "analysis_output/scripts/15_cross_theme_meta_analysis.py")
m15 = importlib.util.module_from_spec(spec15)
spec15.loader.exec_module(m15)

DATA_DIR = "analysis_output/data/"
CRC81861 = "/Users/yun/CRC_raw_data/GSE81861/"
CRC_PAN = "/Users/yun/CRC_raw_data/pan_cancer/"
SDL = "/Users/yun/SignatureDL_raw_data/pan_cancer_phase2/"

TISCH2_PATHS = {
    "GSE167297": CRC_PAN + "STAD_GSE167297/",
    "BLCA_GSE130001": CRC_PAN + "BLCA_GSE130001/",
    "SKCM_GSE120575": SDL + "SKCM_GSE120575/",
    "Glioma_GSE131928": SDL + "Glioma_GSE131928/",
    "AML_GSE116256": SDL + "AML_GSE116256/",
    "KIRC_GSE159115": SDL + "KIRC_GSE159115/",
    "OS_GSE162454": SDL + "OS_GSE162454/",
    "UCEC_GSE139555": SDL + "UCEC_GSE139555/",
    "PAAD_GSE154778": SDL + "PAAD_GSE154778/",
    "NSCLC_GSE131907": SDL + "NSCLC_GSE131907/",
    "CESC_GSE168652": SDL + "CESC_GSE168652/",
    "HNSC_GSE139324": SDL + "HNSC_GSE139324/",
    "OSCC_GSE172577": SDL + "OSCC_GSE172577/",
    "UVM_GSE139829": SDL + "UVM_GSE139829/",
    "ALL_GSE132509": SDL + "ALL_GSE132509/",
    "OV_GSE154600": SDL + "OV_GSE154600/",
    "PRAD_GSE172301": SDL + "PRAD_GSE172301/",
    "MB_GSE155446": SDL + "MB_GSE155446/",
}


def wilcoxon_attribution(expr_df: pd.DataFrame, celltype: pd.Series, genes) -> dict:
    n_total = len(celltype)
    cts = celltype.unique()
    ct_values = celltype.values
    out = {}
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
            sigma_sq = (n1 * n2 / 12.0) * ((n_total + 1) - tie_term / (n_total * (n_total - 1))) if n_total > 1 else 0
            sigma = np.sqrt(sigma_sq) if sigma_sq > 0 else 0
            z = (U1 - mu) / sigma if sigma > 0 else 0.0
            if z > best_z:
                best_z, best_ct = z, ct
        out[g] = best_ct
    return out


def load_gse81861(genes):
    def load(fp):
        df = pd.read_csv(fp, index_col=0)
        df.index = [str(i).split("_")[1] if len(str(i).split("_")) >= 3 else str(i) for i in df.index]
        return df
    nm = load(CRC81861 + "GSE81861_CRC_NM_all_cells_FPKM.csv")
    tm = load(CRC81861 + "GSE81861_CRC_tumor_all_cells_FPKM.csv")
    nm = nm.groupby(nm.index).mean(); tm = tm.groupby(tm.index).mean()
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


def load_brca(genes):
    data_dir = CRC_PAN + "BRCA_GSE176078/Wu_etal_2021_BRCA_scRNASeq/"
    gene_list = pd.read_csv(data_dir + "count_matrix_genes.tsv", header=None, names=["symbol"])["symbol"].tolist()
    barcodes = pd.read_csv(data_dir + "count_matrix_barcodes.tsv", header=None, names=["barcode"])["barcode"].tolist()
    meta = pd.read_csv(data_dir + "metadata.csv", index_col=0)
    mat = mmread(data_dir + "count_matrix_sparse.mtx").tocsr()
    name_to_idx = {}
    for i, n in enumerate(gene_list):
        name_to_idx.setdefault(n, []).append(i)
    present = [g for g in genes if g in name_to_idx]
    total_counts = np.asarray(mat.sum(axis=0)).flatten(); total_counts[total_counts == 0] = 1
    gene_expr = {g: np.asarray(mat[name_to_idx[g], :].sum(axis=0)).flatten() for g in present}
    expr_df = pd.DataFrame(gene_expr, index=barcodes)
    lognorm = np.log1p(expr_df.div(total_counts, axis=0) * 1e4)
    common = lognorm.index.intersection(meta.index)
    lognorm = lognorm.loc[common]
    celltype = pd.Series(meta.loc[common, "celltype_major"].values, index=common)
    return lognorm, celltype


def load_hcc(genes):
    data_dir = CRC_PAN + "HCC_GSE151530/"
    gene_tab = pd.read_csv(data_dir + "GSE151530_genes.tsv.gz", sep="\t", header=None, names=["ensembl", "symbol"])
    barcodes = pd.read_csv(data_dir + "GSE151530_barcodes.tsv.gz", header=None, names=["barcode"])["barcode"].tolist()
    info = pd.read_csv(data_dir + "info.txt", sep="\t").set_index("Cell")
    mat = mmread(data_dir + "GSE151530_matrix.mtx.gz").tocsr()
    name_to_idx = {}
    for i, n in enumerate(gene_tab["symbol"]):
        name_to_idx.setdefault(n, []).append(i)
    present = [g for g in genes if g in name_to_idx]
    gene_expr = {g: np.asarray(mat[name_to_idx[g], :].sum(axis=0)).flatten() for g in present}
    expr_df = pd.DataFrame(gene_expr, index=barcodes)
    total_counts = np.asarray(mat.sum(axis=0)).flatten(); total_counts[total_counts == 0] = 1
    lognorm = np.log1p(expr_df.div(total_counts, axis=0) * 1e4)
    common = lognorm.index.intersection(info.index)
    lognorm = lognorm.loc[common]
    celltype = pd.Series(info.loc[common, "Type"].values, index=common)
    keep = celltype != "unclassified"
    return lognorm.loc[keep], celltype.loc[keep]


def load_tisch2(data_dir, genes):
    f = h5py.File(data_dir + "expression.h5", "r")
    n_genes, n_cells = f["matrix/shape"][:]
    names = np.array([x.decode() for x in f["matrix/features/name"][:]])
    barcodes = np.array([x.decode() for x in f["matrix/barcodes"][:]])
    name_to_idx = {}
    for i, n in enumerate(names):
        name_to_idx.setdefault(n, []).append(i)
    resolved = {}
    for g in genes:
        alt = g.replace(".", "-")
        if g in name_to_idx: resolved[g] = g
        elif alt in name_to_idx: resolved[g] = alt
    data = f["matrix/data"][:]; indices = f["matrix/indices"][:]; indptr = f["matrix/indptr"][:]
    mat = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).tocsr()
    gene_expr = {}
    for g, resolved_name in resolved.items():
        idxs = name_to_idx[resolved_name]
        gene_expr[g] = np.asarray(mat[idxs, :].max(axis=0).todense()).flatten()
    expr_df = pd.DataFrame(gene_expr, index=barcodes)
    meta = pd.read_csv(data_dir + "meta.tsv", sep="\t").set_index("Cell")
    common = expr_df.index.intersection(meta.index)
    expr_df = expr_df.loc[common]
    celltype = pd.Series(meta.loc[common, "Celltype (major-lineage)"].values, index=common)
    return expr_df, celltype


def main():
    zscore_df = pd.read_csv(DATA_DIR + "phase2_pilot_all_combined.csv")
    zscore_df["dataset"] = zscore_df["model"].map(m15.MODEL_TO_DATASET)
    missing_ds = zscore_df[zscore_df["dataset"].isna()]["model"].unique()
    if len(missing_ds):
        print(f"WARNING: {len(missing_ds)} models not in MODEL_TO_DATASET, skipped: {list(missing_ds)}")
    zscore_df = zscore_df.dropna(subset=["dataset"])

    all_wilcoxon_rows = []
    for dataset, sub in zscore_df.groupby("dataset"):
        genes_needed = sorted(sub["gene"].unique())
        n_models = sub["model"].nunique()
        print(f"[{dataset}] {n_models} models, {len(genes_needed)} unique genes ...", flush=True)
        try:
            if dataset == "GSE81861":
                expr, ct = load_gse81861(genes_needed)
            elif dataset == "GSE176078":
                expr, ct = load_brca(genes_needed)
            elif dataset == "GSE151530":
                expr, ct = load_hcc(genes_needed)
            elif dataset in TISCH2_PATHS:
                expr, ct = load_tisch2(TISCH2_PATHS[dataset], genes_needed)
            else:
                print(f"  !! no loader for dataset {dataset}, skipping {n_models} models")
                continue
        except Exception as e:
            print(f"  !! failed to load {dataset}: {e}, skipping {n_models} models")
            continue

        gene_to_ct = wilcoxon_attribution(expr, ct, genes_needed)
        print(f"  -> attributed {len(gene_to_ct)}/{len(genes_needed)} genes")
        for _, row in sub.iterrows():
            ct_call = gene_to_ct.get(row["gene"])
            if ct_call is not None:
                all_wilcoxon_rows.append({"model": row["model"], "gene": row["gene"], "top_celltype": ct_call})

    wilcoxon_df = pd.DataFrame(all_wilcoxon_rows)
    wilcoxon_df.to_csv(DATA_DIR + "phase2_pilot_all_combined_wilcoxon.csv", index=False)
    print(f"\nWilcoxon attribution covers {wilcoxon_df['model'].nunique()}/108 models, "
          f"{len(wilcoxon_df)}/{len(zscore_df)} gene-model rows")

    print("\n" + "=" * 70)
    print("z-score主分析(现有论文结论)")
    print("=" * 70)
    z_res = m15.analyze(zscore_df, exclude_rare=False)
    z_sig = z_res[z_res.significant_fdr05]
    z_sig_genuine = z_sig[~z_sig.tautological_pool]
    z_h1a = z_sig_genuine[z_sig_genuine.top_category.isin(["Stromal_Vascular", "Immune_Hematopoietic"])]
    z_h1b = z_sig_genuine[z_sig_genuine.top_category == "Malignant"]
    print(f"z-score: {len(z_res)}个模型, 显著{len(z_sig)}个({len(z_sig)/len(z_res)*100:.1f}%), "
          f"H1a={len(z_h1a)}, H1b={len(z_h1b)}")

    print("\n" + "=" * 70)
    print("Wilcoxon独立算法重新跑全部模型")
    print("=" * 70)
    w_res = m15.analyze(wilcoxon_df, exclude_rare=False)
    w_sig = w_res[w_res.significant_fdr05]
    w_sig_genuine = w_sig[~w_sig.tautological_pool]
    w_h1a = w_sig_genuine[w_sig_genuine.top_category.isin(["Stromal_Vascular", "Immune_Hematopoietic"])]
    w_h1b = w_sig_genuine[w_sig_genuine.top_category == "Malignant"]
    print(f"Wilcoxon: {len(w_res)}个模型, 显著{len(w_sig)}个({len(w_sig)/len(w_res)*100:.1f}%), "
          f"H1a={len(w_h1a)}, H1b={len(w_h1b)}")

    w_res.to_csv(DATA_DIR + "phase3_meta_analysis_wilcoxon.csv", index=False)

    # 模型级别一致性:两种算法对"该模型是否显著"的判断是否一致
    merged = z_res[["model", "significant_fdr05", "top_category"]].merge(
        w_res[["model", "significant_fdr05", "top_category"]], on="model", suffixes=("_z", "_w"), how="outer")
    merged["sig_agree"] = merged["significant_fdr05_z"] == merged["significant_fdr05_w"]
    merged["category_agree"] = merged["top_category_z"] == merged["top_category_w"]
    print(f"\n模型级别'是否显著'判断一致: {merged['sig_agree'].sum()}/{len(merged)} "
          f"({merged['sig_agree'].sum()/len(merged)*100:.1f}%)")
    print(f"模型级别众数细胞类型大类判断一致: {merged['category_agree'].sum()}/{len(merged)} "
          f"({merged['category_agree'].sum()/len(merged)*100:.1f}%)")
    disagree_sig = merged[~merged.sig_agree]
    print(f"\n判断'是否显著'不一致的模型({len(disagree_sig)}个):")
    print(disagree_sig.to_string(index=False))
    merged.to_csv(DATA_DIR + "wilcoxon_vs_zscore_full_model_comparison.csv", index=False)


if __name__ == "__main__":
    main()
