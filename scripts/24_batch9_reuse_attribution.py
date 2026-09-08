"""
阶段2扩大样本第六批(batch9,来自新主题词表扩检索出的候选池,命中率62.9%的富矿):
16个命中模型癌种和现成数据集(BLCA/STAD/UCEC/CRC/HCC/KIRC/NSCLC/Glioma/BRCA)完全重合,
零新增数据成本。另有6个命中(Wilms瘤、髓母细胞瘤、葡萄膜黑色素瘤×2、口腔鳞癌×2)是全新
癌种,本地无参考数据集,本轮不处理。
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
    "BLCA_GSE130001": {
        "dir": CRC_PAN + "BLCA_GSE130001/",
        "models": {
            "PMC7880321_immune_BLCA": ["AHNAK", "CALR", "PLA2G2A", "OAS1", "CDH1", "PDGFRA",
                                         "SEMA3F", "RAC3", "IRF5", "CARD11"],
            "PMC12891303_CAF_BLCA": ["SERPINF1", "DDR2", "SLIT2", "HSPG2", "ECM1", "RECK"],
        },
    },
    "STAD_GSE167297": {
        "dir": CRC_PAN + "STAD_GSE167297/",
        "models": {
            "PMC11634757_lactylation_STAD": ["COL4A1", "SLC16A7", "IRAK1"],
        },
    },
    "UCEC_GSE139555": {
        "dir": SDL + "UCEC_GSE139555/",
        "models": {
            "PMC9578220_chromatin_immune_UCEC": ["BTNL9", "CD40LG", "CD47", "HLA-DMB", "HLA-DRB5",
                                                   "HLA-G", "TNFRSF14", "TNFRSF18", "TNFRSF4"],
            "PMC12891315_CAF_UCEC": ["ANLN", "ASPM", "CENPI", "CKAP2L", "DNA2", "KIF14", "LTB",
                                       "POLQ", "TTK", "ZNF695", "HAPLN1", "CDK16", "CIT", "FOXL2",
                                       "ISLR", "LRRN4CL", "RARRES2", "SERPING1", "TICRR"],
        },
    },
    "KIRC_GSE159115": {
        "dir": SDL + "KIRC_GSE159115/",
        "models": {
            "PMC9614380_aminoacid_KIRC": ["ACADSB", "BCKDHB", "CBS", "CSAD", "HAO1", "HIBCH",
                                            "HOGA1", "IL4I1", "IYD", "NNMT", "PSAT1", "PYCR1",
                                            "RIMKLA", "RPL13", "RPL22L1", "RPL36A", "SLC5A5", "UROC1"],
        },
    },
    "NSCLC_GSE131907": {
        "dir": SDL + "NSCLC_GSE131907/",
        "models": {
            "PMC11310873_lncRNA_LUAD": ["AC092168.2", "LINC01352", "LINC00968", "AC024075.1",
                                          "AC005070.3", "AL133445.2", "AC005856.1"],
            "PMC10724919_ICD_lncRNA_LUSC": ["MIR22HG", "LINC02345", "AC137932.2", "AP001189.1",
                                              "AC007823.1", "AC087521.1", "AP001189.3", "LRRK2-DT",
                                              "AC008972.2", "LINC02471", "AC009570.2", "LNCOG"],
            "PMC9843944_autophagy_LUAD": ["APOL1", "ARSA", "ATG10", "ATG12", "ATG4A", "BCL2L1",
                                            "BNIP3L", "CAPNS1", "CX3CL1", "DAPK2", "EEF2K", "EIF2S1",
                                            "EIF4G1", "GNAI3", "KLHL24", "MAPK8IP1", "MBTPS2", "NLRC4",
                                            "PRKCD", "RELA", "SIRT2", "SPHK1", "ST13"],
            "PMC8379743_m6A_LUAD": ["CLEC3B", "TENM3", "IGF2BP1", "E2F7", "ANLN", "ANKRD18B", "FBN2"],
        },
    },
    "Glioma_GSE131928": {
        "dir": SDL + "Glioma_GSE131928/",
        "models": {
            "PMC9618960_ICD_LGG": ["IL17RA", "IL1R1", "EIF2AK3", "CD4", "PRF1", "CXCR3", "CD8A",
                                     "BAX", "PDIA3", "CASP8", "MYD88", "CASP1"],
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
        "PMC10932653_chromatin_HCC": ["BMI1", "CBX2", "MRGBP"],
        "PMC10662873_anoikis_HCC": ["BSG", "PLK1", "SPP1", "PBK", "NQO1"],
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

    genes = ["MMP3", "HSPA8", "PTGIS", "CPT2", "GPI", "CDKN2A", "SPP1", "GRP", "APOE", "CHGA",
             "CAV1", "GSTM1", "VEGFA", "LAMA2", "TIMP1", "AGRN", "HSPA1A", "SNAP25"]
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
    sub = pd.DataFrame({"model": "PMC10531437_lipid_ERstress_CRC", "gene": genes_present,
                         "top_celltype": top_ct[genes_present].values,
                         "zscore": z.loc[genes_present].max(axis=1).round(2).values})
    print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
    print("分布:", dict(Counter(sub.top_celltype)))
    return sub


def run_brca():
    print("\n########## GSE176078 (BRCA, 追加2模型) ##########")
    data_dir = CRC_PAN + "BRCA_GSE176078/Wu_etal_2021_BRCA_scRNASeq/"
    models = {
        "PMC11085026_ESR_BRCA": ["ELOVL2", "IFNG", "IGF2BP1", "MAP2K6", "MZB1", "PCSK6", "PCSK9", "POP1"],
        "PMC10791883_CAF_BRCA": ["C1S", "CCDC8", "COL12A1", "CTSO", "IGFBP6", "MFAP4", "NT5E",
                                   "OSMR", "PDLIM4", "RUNX1", "SAV1", "SDC1", "SGCE", "TLN2"],
    }
    gene_list = pd.read_csv(data_dir + "count_matrix_genes.tsv", header=None, names=["symbol"])["symbol"].tolist()
    barcodes = pd.read_csv(data_dir + "count_matrix_barcodes.tsv", header=None, names=["barcode"])["barcode"].tolist()
    meta = pd.read_csv(data_dir + "metadata.csv", index_col=0)
    mat = mmread(data_dir + "count_matrix_sparse.mtx").tocsr()
    name_to_idx = {}
    for i, n in enumerate(gene_list):
        name_to_idx.setdefault(n, []).append(i)
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
    common = lognorm.index.intersection(meta.index)
    lognorm = lognorm.loc[common]
    lognorm["celltype"] = meta.loc[common, "celltype_major"].values
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


def main():
    all_dfs = []
    for tag, cfg in TISCH2_JOBS.items():
        all_dfs.append(run_tisch2(tag, cfg))
    all_dfs.append(run_hcc())
    all_dfs.append(run_crc())
    all_dfs.append(run_brca())
    out = pd.concat(all_dfs, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_batch9_reuse_attribution.csv", index=False)
    print(f"\n写入 {OUT_DIR}phase2_batch9_reuse_attribution.csv, 总模型数: {out['model'].nunique()}, 总基因数: {len(out)}")


if __name__ == "__main__":
    main()
