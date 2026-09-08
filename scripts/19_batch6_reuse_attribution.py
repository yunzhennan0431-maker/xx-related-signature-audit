"""
阶段2扩大样本第三批(batch6):新一轮批量复核命中里,12个模型的癌种和现成数据集
(STAD/SKCM/HCC/Glioma/CESC/NSCLC/KIRC/CRC)完全重合,另有BLCA为新增数据集
(GSE130001,h5+meta.tsv格式已在本地,零下载成本),一并跑归因。
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

# ---- TISCH2格式数据集(h5+meta.tsv) ----
TISCH2_JOBS = {
    "STAD_GSE167297": {
        "dir": CRC_PAN + "STAD_GSE167297/",
        "models": {
            "PMC9767988_pyroptosis_lncRNA_STAD": ["NCRNA00094", "TUG1", "LOC541471", "AC058791.2",
                                                    "JPX", "GNAS-AS1", "MGC12916", "RP11-834C11.4"],
            "PMC9614251_glycolysis_lncRNA_STAD": ["AL353804.1", "AC010719.1", "TNFRSF10A-AS1",
                                                    "AC005586.1", "AC009948.1", "AL355574.1", "AL161785.1"],
        },
    },
    "SKCM_GSE120575": {
        "dir": SDL + "SKCM_GSE120575/",
        "models": {
            "PMC11178724_cuproptosis_methylation_SKCM": ["AOC1", "ACO1", "DPYD", "ABCB8", "CIAO1"],
        },
    },
    "Glioma_GSE131928": {
        "dir": SDL + "Glioma_GSE131928/",
        "models": {
            "PMC9342864_autophagy_Glioma": ["DIRAS3", "CFLAR", "BAX", "TP53", "GRID2", "BIRC5",
                                             "MAPK9", "PTK6", "MYC"],
        },
    },
    "CESC_GSE168652": {
        "dir": SDL + "CESC_GSE168652/",
        "models": {
            "PMC9669765_autophagy_lncRNA_CESC": ["MIR9-3HG", "SMURF2P1", "AC005332.4"],
        },
    },
    "NSCLC_GSE131907": {
        "dir": SDL + "NSCLC_GSE131907/",
        "models": {
            "PMC9561419_glycolysis_immune_LUSC": ["PYGB", "GDF2", "SERPIND1", "MDH1", "TSLP"],
            "PMC10414028_hypoxia_LUAD": ["HK1", "PDK3", "PFKL", "SLC2A1", "STC1", "XPNPEP1"],
            "PMC9550247_autophagy_LUAD": ["EIF2AK3", "ITGB1"],
        },
    },
    "KIRC_GSE159115": {
        "dir": SDL + "KIRC_GSE159115/",
        "models": {
            "PMC8403747_glycolysis_lncRNA_KIRC": ["AC026401.3", "AC087741.1", "AC008906.1", "IGFL2-AS1",
                                                    "ADAMTS9-AS2", "SPINT1-AS1", "ATP1A1-AS1"],
        },
    },
    "BLCA_GSE130001": {
        "dir": CRC_PAN + "BLCA_GSE130001/",
        "models": {
            "PMC11953504_disulfidptosis_lncRNA_BLCA": ["LINC-PINT", "AC023825.2", "AC010331.1", "AL590428.1",
                                                          "LSAMP-AS1", "AC009716.1", "AC104785.1", "AC008764.6",
                                                          "LINC01184"],
            # 注:原模型共12变量,另3个(T_cells_CD8/Macrophages_M1/Eosinophils)是CIBERSORT
            # 免疫细胞比例反卷积特征而非基因,不纳入单细胞基因归因
            "PMC9341065_immune_autophagy_BLCA": ["C5AR2", "CD96", "CSF3R", "FBXW10", "FCAR", "GHR",
                                                   "IL10", "MEFV", "OLR1", "PGLYRP3", "RASGRP4", "S100A12"],
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
    return pd.concat(rows, ignore_index=True), dict(meta["Celltype (major-lineage)"].value_counts())


def run_hcc():
    print("\n########## GSE151530 (HCC, 追加2模型) ##########")
    CRC = CRC_PAN + "HCC_GSE151530/"
    genes = pd.read_csv(CRC + "GSE151530_genes.tsv.gz", sep="\t", header=None, names=["ensembl", "symbol"])
    barcodes = pd.read_csv(CRC + "GSE151530_barcodes.tsv.gz", header=None, names=["barcode"])["barcode"].tolist()
    info = pd.read_csv(CRC + "info.txt", sep="\t").set_index("Cell")
    mat = mmread(CRC + "GSE151530_matrix.mtx.gz").tocsr()
    name_to_idx = {}
    for i, n in enumerate(genes["symbol"]):
        name_to_idx.setdefault(n, []).append(i)

    models = {
        "PMC11382408_cuproptosis_HCC": ["SOX4", "NDRG2", "MYC", "TM4SF1", "CYB5A", "IFI27"],
        "PMC8605142_glycolysis_HCC": ["CDCA8", "RAB5IF", "SAP30", "UCK2"],
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

    genes = ["PPP2CB", "PGM2", "PPARGC1A", "ENO3", "PMM2", "P4HA1", "STC2", "CHPF2"]
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
    sub = pd.DataFrame({"model": "PMC8377503_glycolysis_colon", "gene": genes_present,
                         "top_celltype": top_ct[genes_present].values,
                         "zscore": z.loc[genes_present].max(axis=1).round(2).values})
    print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
    print("分布:", dict(Counter(sub.top_celltype)))
    return sub


def main():
    all_dfs = []
    dataset_counts = {}
    for tag, cfg in TISCH2_JOBS.items():
        df, counts = run_tisch2(tag, cfg)
        all_dfs.append(df)
        dataset_counts[tag] = counts
    all_dfs.append(run_hcc())
    all_dfs.append(run_crc())
    out = pd.concat(all_dfs, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_batch6_reuse_attribution.csv", index=False)
    print(f"\n写入 {OUT_DIR}phase2_batch6_reuse_attribution.csv, 总模型数: {out['model'].nunique()}, 总基因数: {len(out)}")
    print("\n各新数据集候选细胞类型及细胞数(供更新DATASET_CELLTYPE_COUNTS):")
    for tag, counts in dataset_counts.items():
        print(tag, ":", counts)


if __name__ == "__main__":
    main()
