"""
batch11(继续扩大候选池复核,命中率25.7%,含1处重要发现:PMC13114778被证实疑似重复使用
PMC10900655的模型输出并误标癌种,见15_cross_theme_meta_analysis.py的MODEL_TO_DATASET注释)
的6个零新增数据成本命中模型归因脚本。另有1个命中(PMC12973827,神经母细胞瘤)本地无参考
单细胞数据,留待后续。
"""
import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sp
from collections import Counter

OUT_DIR = "analysis_output/data/"
SDL = "/Users/yun/SignatureDL_raw_data/pan_cancer_phase2/"
CRC_PAN = "/Users/yun/CRC_raw_data/pan_cancer/"
CRC81861 = "/Users/yun/CRC_raw_data/GSE81861/"

TISCH2_JOBS = {
    "OV_GSE154600": {
        "dir": SDL + "OV_GSE154600/",
        "models": {
            "PMC10807923_anoikis_OV": ["NTRK2", "FN1", "MTOR", "IGF1", "ABHD4", "HMCN1", "CDKN1B",
                                         "KRAS", "MAVS", "RB1", "CRYAB", "LRP1", "LPAR1"],
        },
    },
    "PRAD_GSE172301": {
        "dir": SDL + "PRAD_GSE172301/",
        "models": {
            "PMC10469923_anoikis_PRAD": ["EGF", "MYC", "PLK1", "EZH2", "AFP", "NOX4", "BMP6", "MMP11"],
        },
    },
    "HNSC_GSE139324": {
        "dir": SDL + "HNSC_GSE139324/",
        "models": {
            "PMC10683111_anoikis_lncRNA_HNSCC": ["CYTOR", "Z97200.1", "RAB11B-AS1", "LINC02084",
                                                    "EMSLR", "LINC01503", "AC015878.1"],
        },
    },
}


def run_tisch2(tag, cfg):
    data_dir = cfg["dir"]
    models = cfg["models"]
    print(f"\n########## {tag} ##########")
    f = h5py.File(data_dir + "expression.h5", "r")
    n_genes, n_cells = f["matrix/shape"][:]
    names = np.array([x.decode() for x in f["matrix/features/name"][:]])
    barcodes = np.array([x.decode() for x in f["matrix/barcodes"][:]])
    name_to_idx = {}
    for i, n in enumerate(names):
        name_to_idx.setdefault(n, []).append(i)
    all_raw_genes = sorted(set(g for genes in models.values() for g in genes))
    resolved = {}
    for g in all_raw_genes:
        alt = g.replace(".", "-")
        alt2 = g.replace("-", ".")
        if g in name_to_idx: resolved[g] = g
        elif alt in name_to_idx: resolved[g] = alt
        elif alt2 in name_to_idx: resolved[g] = alt2
    present = list(resolved.values())
    print("present:", len(present), "/", len(all_raw_genes), "missing:", [g for g in all_raw_genes if g not in resolved])

    data = f["matrix/data"][:]; indices = f["matrix/indices"][:]; indptr = f["matrix/indptr"][:]
    mat = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).tocsr()
    gene_expr = {}
    for g in present:
        idxs = name_to_idx[g]
        gene_expr[g] = np.asarray(mat[idxs, :].max(axis=0).todense()).flatten()
    expr_df = pd.DataFrame(gene_expr, index=barcodes)
    meta = pd.read_csv(data_dir + "meta.tsv", sep="\t").set_index("Cell")
    ct_col = "Celltype (major-lineage)" if "Celltype (major-lineage)" in meta.columns else [
        c for c in meta.columns if "major" in c.lower() or "celltype" in c.lower()
    ][0]
    common = expr_df.index.intersection(meta.index)
    expr_df = expr_df.loc[common]
    expr_df["celltype"] = meta.loc[common, ct_col].values
    mean_by_ct = expr_df.groupby("celltype")[present].mean().T
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)

    rows = []
    for model_name, raw_genes in models.items():
        mp = [resolved[g] for g in raw_genes if g in resolved]
        sub = pd.DataFrame({"model": model_name, "gene": mp, "top_celltype": top_ct[mp].values,
                             "zscore": z.loc[mp].max(axis=1).round(2).values})
        rows.append(sub)
        print(f"\n=== {model_name} ({len(mp)}/{len(raw_genes)}) ===")
        print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
        print("分布:", dict(Counter(sub.top_celltype)))
        print("候选细胞类型(K):", len(mean_by_ct.columns))
    return pd.concat(rows, ignore_index=True)


def run_crc():
    print("\n########## GSE81861 (CRC, 追加2模型) ##########")

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
    cols = combined.columns
    cell_type = [c.split("__")[1] if "__" in c else "NA" for c in cols]
    meta = pd.DataFrame({"cell": cols, "celltype": cell_type})
    meta = meta[meta.celltype != "NA"]
    combined = combined[meta.cell.values]

    models = {
        "PMC9465161_m7G_lncRNA_COAD": ["ITFG1-AS1", "ATP2B1-AS1", "LINC02257", "SEPTIN7-DT",
                                         "LINC02593", "NSMCE1-DT", "LINC01011", "PRKAR1B-AS2",
                                         "ALMS1-IT1", "ALKBH3-AS1", "LENG8-AS1", "NDUFB2-AS1",
                                         "LINC01909", "LINC02428"],
        "PMC10932661_m7G_lncRNA_mucinousCRC": ["AC090152.1", "AC254629.1", "LINC01133", "LINC01134",
                                                  "MAN1B1-DT", "MHENCR", "SMIM2-AS1", "XACT"],
    }
    all_raw_genes = sorted(set(g for genes_ in models.values() for g in genes_))
    present_all = [g for g in all_raw_genes if g in combined.index]
    print("present:", len(present_all), "/", len(all_raw_genes), "missing:", [g for g in all_raw_genes if g not in combined.index])
    logexpr = np.log2(combined.loc[present_all] + 1)
    result = {}
    for ct in meta.celltype.unique():
        cells = meta[meta.celltype == ct].cell.values
        result[ct] = logexpr[cells].mean(axis=1)
    mean_by_ct = pd.DataFrame(result)
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)

    rows = []
    for model_name, raw_genes in models.items():
        mp = [g for g in raw_genes if g in present_all]
        sub = pd.DataFrame({"model": model_name, "gene": mp, "top_celltype": top_ct[mp].values,
                             "zscore": z.loc[mp].max(axis=1).round(2).values})
        rows.append(sub)
        print(f"\n=== {model_name} ({len(mp)}/{len(raw_genes)}) ===")
        print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
        print("分布:", dict(Counter(sub.top_celltype)))
    return pd.concat(rows, ignore_index=True)


def run_brca():
    print("\n########## GSE176078 (BRCA, 追加1模型) ##########")
    data_dir = CRC_PAN + "BRCA_GSE176078/Wu_etal_2021_BRCA_scRNASeq/"
    genes = ["YAP1", "PIK3R1", "BAK1", "PHLDA2", "EDA2R", "LAMB3", "CD24", "SLC2A1", "CDC25C", "SLC39A6"]
    gene_list = pd.read_csv(data_dir + "count_matrix_genes.tsv", header=None, names=["symbol"])["symbol"].tolist()
    barcodes = pd.read_csv(data_dir + "count_matrix_barcodes.tsv", header=None, names=["barcode"])["barcode"].tolist()
    meta = pd.read_csv(data_dir + "metadata.csv", index_col=0)
    from scipy.io import mmread
    mat = mmread(data_dir + "count_matrix_sparse.mtx").tocsr()
    name_to_idx = {}
    for i, n in enumerate(gene_list):
        name_to_idx.setdefault(n, []).append(i)
    present = [g for g in genes if g in name_to_idx]
    print("present:", len(present), "/", len(genes), "missing:", [g for g in genes if g not in name_to_idx])

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
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)
    sub = pd.DataFrame({"model": "PMC11411761_anoikis_BRCA", "gene": present,
                         "top_celltype": top_ct[present].values,
                         "zscore": z.loc[present].max(axis=1).round(2).values})
    print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
    print("分布:", dict(Counter(sub.top_celltype)))
    return sub


def run_stad():
    print("\n########## GSE167297 (STAD, 追加1模型) ##########")
    data_dir = CRC_PAN + "STAD_GSE167297/"
    genes = ["MAGI2-AS3", "LINC00106", "AC145285.6", "AL590705.3", "AC007405.3"]
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
    present = list(resolved.values())
    print("present:", len(present), "/", len(genes), "missing:", [g for g in genes if g not in resolved])
    data = f["matrix/data"][:]; indices = f["matrix/indices"][:]; indptr = f["matrix/indptr"][:]
    mat = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).tocsr()
    gene_expr = {}
    for g in present:
        idxs = name_to_idx[g]
        gene_expr[g] = np.asarray(mat[idxs, :].max(axis=0).todense()).flatten()
    expr_df = pd.DataFrame(gene_expr, index=barcodes)
    meta = pd.read_csv(data_dir + "meta.tsv", sep="\t").set_index("Cell")
    ct_col = "Celltype (major-lineage)"
    common = expr_df.index.intersection(meta.index)
    expr_df = expr_df.loc[common]
    expr_df["celltype"] = meta.loc[common, ct_col].values
    mean_by_ct = expr_df.groupby("celltype")[present].mean().T
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)
    sub = pd.DataFrame({"model": "PMC10339815_DDR_lncRNA_STAD", "gene": present,
                         "top_celltype": top_ct[present].values,
                         "zscore": z.loc[present].max(axis=1).round(2).values})
    print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
    print("分布:", dict(Counter(sub.top_celltype)))
    return sub


def main():
    all_dfs = []
    for tag, cfg in TISCH2_JOBS.items():
        all_dfs.append(run_tisch2(tag, cfg))
    all_dfs.append(run_crc())
    all_dfs.append(run_brca())
    all_dfs.append(run_stad())
    out = pd.concat(all_dfs, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_batch11_reuse_attribution.csv", index=False)
    print(f"\n写入 {OUT_DIR}phase2_batch11_reuse_attribution.csv, 总模型数: {out['model'].nunique()}, 总基因数: {len(out)}")


if __name__ == "__main__":
    main()
