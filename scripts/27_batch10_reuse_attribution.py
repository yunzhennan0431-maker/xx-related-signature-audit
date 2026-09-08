"""
阶段2扩大样本第七批(batch10,继续深挖主题词表扩充后的候选池,命中率48.6%):
15个命中模型癌种和现成数据集(Glioma/HCC/SKCM/OS/BLCA/UCEC/NSCLC/CRC/STAD/BRCA)完全重合,
零新增数据成本。另有2个命中(头颈鳞癌NK-TAM互作、急性淋巴细胞白血病FBXO8)是全新癌种,
本地无参考数据集,本轮不处理。
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
    "Glioma_GSE131928": {
        "dir": SDL + "Glioma_GSE131928/",
        "models": {
            "PMC12354743_NETs_Glioma": ["MAPK1", "P2RX1", "PARVB", "STAT3"],
        },
    },
    "SKCM_GSE120575": {
        "dir": SDL + "SKCM_GSE120575/",
        "models": {
            "PMC10980369_ERstress_SKCM": ["ZBP1", "DIABLO", "GNLY", "FASLG", "AURKA", "TNFRSF21", "CD40LG"],
            "PMC9361349_DNArepair_SKCM": ["TYMS", "SNAPC5", "CMPK2", "PDE4B", "HCLS1", "NME1", "POLR2A",
                                           "COX17", "LIG1", "POLE4", "GTF2H1", "AK1"],
            "PMC12179157_PANoptosis_SKCM": ["CCL4", "CTSS", "ADAMDEC1", "CMAHP", "CD69", "PIM2",
                                              "CRIP1", "LSP1", "BCL11B", "CCR7"],
        },
    },
    "OS_GSE162454": {
        "dir": SDL + "OS_GSE162454/",
        "models": {
            "PMC10076677_anoikis_lncRNA_OS": ["AC079612.1", "MEF2C-AS1", "SNHG6", "TBX2-AS1"],
        },
    },
    "BLCA_GSE130001": {
        "dir": CRC_PAN + "BLCA_GSE130001/",
        "models": {
            "PMC13355270_mitophagy_BLCA": ["NTN4", "AKR1B1", "FBXW7", "MYH10", "IGF2BP3"],
        },
    },
    "UCEC_GSE139555": {
        "dir": SDL + "UCEC_GSE139555/",
        "models": {
            "PMC12055729_lactylation_UCEC": ["ZFHX4", "SCGB2A1", "IGSF1"],
        },
    },
    "NSCLC_GSE131907": {
        "dir": SDL + "NSCLC_GSE131907/",
        "models": {
            "PMC9870742_anoikis_LUAD": ["CDH2", "PARP1", "LAMB3", "ANGPTL4", "PPARG", "PAK1", "CAV1",
                                          "PLAT", "FGF2", "ITGA8", "ITGB4", "HMGA1", "SLC2A1", "PLK1",
                                          "TIMP1", "PBK"],
        },
    },
    "STAD_GSE167297": {
        "dir": CRC_PAN + "STAD_GSE167297/",
        "models": {
            "PMC8291813_cGASSTING_STAD": ["IFNB1", "IFNA4", "IL6", "NFKB2", "TRIM25"],
        },
    },
    "KIRC_GSE159115": {
        "dir": SDL + "KIRC_GSE159115/",
        "models": {
            "PMC9523360_telomere_KIRC": ["ISG15", "RFC2", "TRIM15", "NEK6", "PRKCQ", "ATP1A1",
                                          "ELOVL3", "TUBB2B", "PLCL1", "NR1H3"],
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
    CRC = CRC_PAN + "HCC_GSE151530/"
    genes = pd.read_csv(CRC + "GSE151530_genes.tsv.gz", sep="\t", header=None, names=["ensembl", "symbol"])
    barcodes = pd.read_csv(CRC + "GSE151530_barcodes.tsv.gz", header=None, names=["barcode"])["barcode"].tolist()
    info = pd.read_csv(CRC + "info.txt", sep="\t").set_index("Cell")
    mat = mmread(CRC + "GSE151530_matrix.mtx.gz").tocsr()
    name_to_idx = {}
    for i, n in enumerate(genes["symbol"]):
        name_to_idx.setdefault(n, []).append(i)

    models = {
        "PMC9468367_ICD_lncRNA_HCC": ["AL049840.5", "AL096678.1", "AC099066.2", "MKLN1-AS", "AC091057.3",
                                        "LINC01224", "AL158163.1", "AC034229.4", "DUXAP8", "GABPB1-AS1"],
        "PMC9720162_ICD_lncRNA_HCC": ["TMEM220-AS1", "LINC02362", "LINC01554", "LINC02499"],
        "PMC10616163_chromatin_HCC": ["PPARGC1A", "DUSP1", "APOBEC3A", "AIRE", "HDAC11", "HMGB2", "APOBEC3B"],
        "PMC11866598_telomere_HCC": ["ASF1A", "CDCA8", "HMMR", "IPO13", "MT3", "MYCN", "PPM1G", "RAD54B",
                                       "RRAGC", "RTN3", "RUVBL1", "SAP30", "SLC7A11", "SMG5", "TALDO1",
                                       "TCOF1", "TRAPPC4", "UAP1L1"],
    }
    all_raw_genes = sorted(set(g for genes_ in models.values() for g in genes_))
    resolved = {}
    for g in all_raw_genes:
        alt = g.replace(".", "-")
        if g in name_to_idx: resolved[g] = g
        elif alt in name_to_idx: resolved[g] = alt
    present = list(resolved.values())
    print("present:", len(present), "/", len(all_raw_genes), "missing:", [g for g in all_raw_genes if g not in resolved])

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
        mp = [resolved[g] for g in raw_genes if g in resolved]
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

    genes = ["AMBRA1", "ATG14", "MAP1LC3A", "MAP1LC3B", "OPTN", "VDAC1", "ATG5", "CSNK2A2", "MFN1", "TOMM22"]
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
    sub = pd.DataFrame({"model": "PMC9636133_mitophagy_CRC", "gene": genes_present,
                         "top_celltype": top_ct[genes_present].values,
                         "zscore": z.loc[genes_present].max(axis=1).round(2).values})
    print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
    print("分布:", dict(Counter(sub.top_celltype)))
    return sub


def run_brca():
    print("\n########## GSE176078 (BRCA, 追加1模型) ##########")
    data_dir = CRC_PAN + "BRCA_GSE176078/Wu_etal_2021_BRCA_scRNASeq/"
    genes = ["HSPA1B", "LGALS4", "HSPB1", "ARPC1B", "ACTG1", "LEPROTL1", "HMGN2", "SFPQ", "RBM3",
             "SLC2A3", "LY6E", "ARL4C", "UBE2I", "AP2M1", "TERF2IP", "H2AFY", "RAB5C", "LITAF",
             "MED10", "RBM17", "METTL9", "NR1H2", "PTTG1"]
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
    sub = pd.DataFrame({"model": "PMC13114778_lactylation_BRCA", "gene": present,
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
    all_dfs.append(run_crc())
    all_dfs.append(run_brca())
    out = pd.concat(all_dfs, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_batch10_reuse_attribution.csv", index=False)
    print(f"\n写入 {OUT_DIR}phase2_batch10_reuse_attribution.csv, 总模型数: {out['model'].nunique()}, 总基因数: {len(out)}")


if __name__ == "__main__":
    main()
