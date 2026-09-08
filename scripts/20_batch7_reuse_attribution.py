"""
阶段2扩大样本第四批(batch7):437篇候选池里score降到3-4区间后的下一批复核,
14个命中模型癌种和现成数据集(OS/KIRC/PAAD/Glioma/STAD/BRCA/CRC/HCC/UCEC/NSCLC)
完全重合,零新增数据成本。
"""
import numpy as np
import pandas as pd
from scipy.io import mmread
import h5py
import scipy.sparse as sp
from collections import Counter

OUT_DIR = "analysis_output/data/"
SDL = "/Users/yun/SignatureDL_raw_data/pan_cancer_phase2/"
CRC_PAN = "/Users/yun/CRC_raw_data/pan_cancer/"
CRC81861 = "/Users/yun/CRC_raw_data/GSE81861/"

TISCH2_JOBS = {
    "OS_GSE162454": {
        "dir": SDL + "OS_GSE162454/",
        "models": {
            "PMC9805482_autophagy_OS": ["CXCR4", "MBTPS2", "MYC"],
        },
    },
    "KIRC_GSE159115": {
        "dir": SDL + "KIRC_GSE159115/",
        "models": {
            "PMC9662042_autophagy_lncRNA_KIRC": ["AL021707.6", "HOTAIRM1", "AC084876.1", "AC010973.2",
                                                   "LINC01507", "AC016773.1", "AC026401.3"],
        },
    },
    "PAAD_GSE154778": {
        "dir": SDL + "PAAD_GSE154778/",
        "models": {
            "PMC9124146_autophagy_lncRNA_PAAD": ["FLVCR1-DT", "AC245041.2", "AC006504.7", "AC125494.2",
                                                   "AC012306.2", "ST20-AS1", "AC036176.1", "LINC01089",
                                                   "AC005696.1", "LINC02257"],
        },
    },
    "Glioma_GSE131928": {
        "dir": SDL + "Glioma_GSE131928/",
        "models": {
            "PMC9097333_autophagy_GBM": ["NDUFB9", "BAK1", "SUPT3H", "GAPDH", "CDKN1B", "CHMP6", "EGFR"],
        },
    },
    "STAD_GSE167297": {
        "dir": CRC_PAN + "STAD_GSE167297/",
        "models": {
            # 注:该模型实际是ESCA+STAD+LIHC+CRC四种消化系统癌种联合建模(18个候选HRlncRNA取
            # 交集后LASSO筛选),不对应单一癌种,这里用STAD作代表性单细胞参考做归因,只是四个
            # 来源癌种之一,解读需要注明这个方法学错配,不能当作和其他单癌种模型同等可比。
            "PMC9092832_hypoxia_lncRNA_digestive_panCancer": ["LUCAT1", "MIR4435-2HG", "LINC01711",
                                                                "AP000695.2", "ADAMTS9-AS2", "AC087521.1"],
            "PMC11913695_disulfidptosis_lncRNA_STAD": ["AC107021.2", "AC129507.1", "FRMD6-AS2"],
            "PMC10046686_pyroptosis_lncRNA_STAD": ["LINC01315", "AL161785.1", "AP003392.1", "AP000695.2",
                                                     "HAGLR", "AL590666.2"],
        },
    },
    "NSCLC_GSE131907": {
        "dir": SDL + "NSCLC_GSE131907/",
        "models": {
            "PMC12210501_ferroptosis_LUAD": ["PPP1R14B", "PLEK2", "RHOV", "C1QTNF6"],
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
        if g in name_to_idx: resolved[g] = g
        elif alt in name_to_idx: resolved[g] = alt
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
    common = expr_df.index.intersection(meta.index)
    expr_df = expr_df.loc[common]
    expr_df["celltype"] = meta.loc[common, "Celltype (major-lineage)"].values
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
    return pd.concat(rows, ignore_index=True)


def run_hcc():
    print("\n########## GSE151530 (HCC, 追加1模型) ##########")
    CRC = CRC_PAN + "HCC_GSE151530/"
    genes = pd.read_csv(CRC + "GSE151530_genes.tsv.gz", sep="\t", header=None, names=["ensembl", "symbol"])
    barcodes = pd.read_csv(CRC + "GSE151530_barcodes.tsv.gz", header=None, names=["barcode"])["barcode"].tolist()
    info = pd.read_csv(CRC + "info.txt", sep="\t").set_index("Cell")
    mat = mmread(CRC + "GSE151530_matrix.mtx.gz").tocsr()
    name_to_idx = {}
    for i, n in enumerate(genes["symbol"]):
        name_to_idx.setdefault(n, []).append(i)

    models = {
        "PMC8100457_lncRNA_HCC": ["AL031985.3", "AL365203.2", "MIR4435-2HG", "AC015908.3"],
    }
    all_raw_genes = sorted(set(g for genes_ in models.values() for g in genes_))
    present = [g for g in all_raw_genes if g in name_to_idx]
    print("present:", len(present), "/", len(all_raw_genes), "missing:", [g for g in all_raw_genes if g not in name_to_idx])

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
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)

    rows = []
    for model_name, raw_genes in models.items():
        mp = [g for g in raw_genes if g in present]
        sub = pd.DataFrame({"model": model_name, "gene": mp, "top_celltype": top_ct[mp].values,
                             "zscore": z.loc[mp].max(axis=1).round(2).values})
        rows.append(sub)
        print(f"\n=== {model_name} ({len(mp)}/{len(raw_genes)}) ===")
        print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
        print("分布:", dict(Counter(sub.top_celltype)))
    return pd.concat(rows, ignore_index=True)


def run_ucec():
    print("\n########## UCEC_GSE139555 (追加1模型) ##########")
    data_dir = SDL + "UCEC_GSE139555/"
    genes = ["AL121906.2", "BOLA3-AS1", "LINC01833", "AC016405.3", "RAB11B-AS1"]
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
    common = expr_df.index.intersection(meta.index)
    expr_df = expr_df.loc[common]
    expr_df["celltype"] = meta.loc[common, "Celltype (major-lineage)"].values
    mean_by_ct = expr_df.groupby("celltype")[present].mean().T
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)
    sub = pd.DataFrame({"model": "PMC7847640_glycolysis_lncRNA_UCEC", "gene": present,
                         "top_celltype": top_ct[present].values,
                         "zscore": z.loc[present].max(axis=1).round(2).values})
    print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
    print("分布:", dict(Counter(sub.top_celltype)))
    return sub


def run_crc():
    print("\n########## GSE81861 (CRC, 追加3模型) ##########")

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
        "PMC8812245_glycolysis_READ": ["ANKZF1", "STC2", "SUCLG2P2", "P4HA1", "GPC1", "PCK1"],
        "PMC11483459_cuproptosis_CRC": ["CDKN2A", "DLAT"],
        "PMC11920466_COAD": ["HEYL", "FSTL3", "CALB2", "HOMER3", "FABP4", "MEGF6", "ADAM8", "TIMP1", "CLCA1", "NOS2"],
    }
    all_raw_genes = sorted(set(g for genes_ in models.values() for g in genes_))
    genes_present = [g for g in all_raw_genes if g in combined.index]
    print("present:", len(genes_present), "/", len(all_raw_genes), "missing:", [g for g in all_raw_genes if g not in combined.index])
    logexpr = np.log2(combined.loc[genes_present] + 1)
    result = {}
    for ct in meta.celltype.unique():
        cells = meta[meta.celltype == ct].cell.values
        result[ct] = logexpr[cells].mean(axis=1)
    mean_by_ct = pd.DataFrame(result)
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)

    rows = []
    for model_name, raw_genes in models.items():
        mp = [g for g in raw_genes if g in genes_present]
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
    genes = ["EIF4EBP1", "IFNG", "TP63"]
    gene_list = pd.read_csv(data_dir + "count_matrix_genes.tsv", header=None, names=["symbol"])["symbol"].tolist()
    barcodes = pd.read_csv(data_dir + "count_matrix_barcodes.tsv", header=None, names=["barcode"])["barcode"].tolist()
    meta = pd.read_csv(data_dir + "metadata.csv", index_col=0)
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
    sub = pd.DataFrame({"model": "PMC9035888_autophagy_BRCA", "gene": present,
                         "top_celltype": top_ct[present].values,
                         "zscore": z.loc[present].max(axis=1).round(2).values})
    print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
    print("分布:", dict(Counter(sub.top_celltype)))
    return sub


def main():
    all_dfs = []
    for tag, cfg in TISCH2_JOBS.items():
        all_dfs.append(run_tisch2(tag, cfg))
    all_dfs.append(run_hcc())
    all_dfs.append(run_ucec())
    all_dfs.append(run_crc())
    all_dfs.append(run_brca())
    out = pd.concat(all_dfs, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_batch7_reuse_attribution.csv", index=False)
    print(f"\n写入 {OUT_DIR}phase2_batch7_reuse_attribution.csv, 总模型数: {out['model'].nunique()}, 总基因数: {len(out)}")


if __name__ == "__main__":
    main()
