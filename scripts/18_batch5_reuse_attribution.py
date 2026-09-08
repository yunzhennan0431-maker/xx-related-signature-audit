"""
阶段2扩大样本第二批:扩大检索命中率76%的一批hit里,有12个模型的癌种和现成数据集
(NSCLC/AML/Glioma/HCC/OS/CRC/SKCM)完全重合,零新增数据成本,直接跑归因。
"""
import numpy as np
import pandas as pd
from scipy.io import mmread
import h5py
import scipy.sparse as sp
from collections import Counter

OUT_DIR = "analysis_output/data/"
SDL = "/Users/yun/SignatureDL_raw_data/pan_cancer_phase2/"
CRC = "/Users/yun/CRC_raw_data/pan_cancer/HCC_GSE151530/"
CRC81861 = "/Users/yun/CRC_raw_data/GSE81861/"

# ---- TISCH2格式数据集(h5) ----
TISCH2_JOBS = {
    "NSCLC_GSE131907": {
        "dir": SDL + "NSCLC_GSE131907/",
        "models": {
            "PMC10085621_hypoxia_NSCLC": ["IGFBP3", "DDIT4", "PHLDA2", "RRAS", "WDR4", "EXO1",
                                           "ECT2", "TYMS", "CDC25C", "PLK1", "FERMT1", "PCSK9"],
            "PMC12268875_disulfidptosis_lncRNA_NSCLC": ["AC096733.2", "ATXN1-AS1", "AF131215.5",
                                                          "AC027288.1", "AL049836.1", "AL606489.1",
                                                          "LINC01711", "AC018645.3"],
        },
    },
    "AML_GSE116256": {
        "dir": SDL + "AML_GSE116256/",
        "models": {
            "PMC11564170_cuproptosis_AML": ["LGR4", "FERP1", "C2orf88", "PGAP1", "ACSM1", "PRSS2",
                                             "IGHD4-17", "SNORD19C", "PSMD6-AS2", "MIR553", "SNRPGP4", "CNN3"],
            "PMC12084325_disulfidptosis_lncRNA_AML": ["AC005076.1", "HDAC4-AS1", "LINC01694",
                                                        "THAP9-AS1", "AP002807.1", "L3MBTL4-AS1"],
        },
    },
    "Glioma_GSE131928": {
        "dir": SDL + "Glioma_GSE131928/",
        "models": {
            "PMC9086515_ferroptosis_Glioma": ["DPP4", "CAPG", "HSPB1", "AURKA", "SESN2", "PGD",
                                               "ARNTL", "LAMP2", "EIF2AK4"],
        },
    },
    "SKCM_GSE120575": {
        "dir": SDL + "SKCM_GSE120575/",
        "models": {
            "PMC9048552_autophagy_lncRNA_SKCM": ["LINC01943", "AC090948.3", "USP30-AS1", "AC068282.1",
                                                   "AC004687.1", "AL133371.2", "AC242842.1", "PCED1B-AS1",
                                                   "HLA-DQB1-AS1", "AC011374.2", "LINC00324", "ITGB2-AS1",
                                                   "AC018553.1", "LINC00520", "DBH-AS1"],
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
    print("\n########## GSE151530 (HCC, 追加4模型) ##########")
    genes = pd.read_csv(CRC + "GSE151530_genes.tsv.gz", sep="\t", header=None, names=["ensembl", "symbol"])
    barcodes = pd.read_csv(CRC + "GSE151530_barcodes.tsv.gz", header=None, names=["barcode"])["barcode"].tolist()
    info = pd.read_csv(CRC + "info.txt", sep="\t").set_index("Cell")
    mat = mmread(CRC + "GSE151530_matrix.mtx.gz").tocsr()
    name_to_idx = {}
    for i, n in enumerate(genes["symbol"]):
        name_to_idx.setdefault(n, []).append(i)

    models = {
        "PMC9524961_hypoxia_HCC": ["ENO1", "SAP30", "STC2"],
        "PMC10294479_autophagy_HCC": ["GNAI3", "FKBP1A", "BIRC5", "SH3GLB1", "HIF1A", "RHEB",
                                       "EIF2S1", "RAB1A", "ATIC", "NPC1", "PRKCD", "ATG4B", "CLN3"],
        "PMC10164321_pyroptosis_HCC": ["BAK1", "GSDME", "NLRP6", "NOD2"],
        "PMC11491388_disulfidptosis_lncRNA_HCC": ["AC009283.1", "TMCC1-AS1", "FOXD2-AS1", "LINC01063", "SLC25A30-AS1"],
    }
    all_raw_genes = sorted(set(g for genes_ in models.values() for g in genes_))
    present = [g for g in all_raw_genes if g in name_to_idx]
    print("present:", len(present), "/", len(all_raw_genes), "missing:", [g for g in all_raw_genes if g not in name_to_idx])

    total_counts = np.asarray(mat.sum(axis=0)).flatten()
    total_counts[total_counts == 0] = 1
    gene_expr = {}
    for g in present:
        idxs = name_to_idx[g]
        gene_expr[g] = np.asarray(mat[idxs, :].sum(axis=0)).flatten()
    expr_df = pd.DataFrame(gene_expr, index=barcodes)
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


def run_crc():
    print("\n########## GSE81861 (CRC, 追加1模型) ##########")

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

    genes = ["VSIG4", "SCG2", "INHBB", "DDC", "CXCL13", "KLK10", "CXCL10", "CCL11"]
    genes_present = [g for g in genes if g in combined.index]
    print("present:", len(genes_present), "/", len(genes), "missing:", [g for g in genes if g not in combined.index])
    logexpr = np.log2(combined.loc[genes_present] + 1)
    result = {}
    for ct in meta.celltype.unique():
        cells = meta[meta.celltype == ct].cell.values
        result[ct] = logexpr[cells].mean(axis=1)
    mean_by_ct = pd.DataFrame(result)
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)
    sub = pd.DataFrame({"model": "PMC11914417_disulfidptosis_CRC", "gene": genes_present,
                         "top_celltype": top_ct[genes_present].values,
                         "zscore": z.loc[genes_present].max(axis=1).round(2).values})
    print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
    print("分布:", dict(Counter(sub.top_celltype)))
    return sub


def main():
    all_dfs = []
    for tag, cfg in TISCH2_JOBS.items():
        all_dfs.append(run_tisch2(tag, cfg))
    all_dfs.append(run_hcc())
    all_dfs.append(run_crc())
    out = pd.concat(all_dfs, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_batch5_reuse_attribution.csv", index=False)
    print(f"\n写入 {OUT_DIR}phase2_batch5_reuse_attribution.csv, 总模型数: {out['model'].nunique()}, 总基因数: {len(out)}")


if __name__ == "__main__":
    main()
