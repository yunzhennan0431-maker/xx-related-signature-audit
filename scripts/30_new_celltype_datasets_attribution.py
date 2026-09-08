"""
弥补"15个新癌种候选一直悬着"这一局限:新下载8个TISCH2数据集,覆盖此前因本地无参考
单细胞数据而搁置的11个模型(另有2个:肾上腺皮质癌/PMC11457769基因数太少(2个,n<6功效不足
且无数据集覆盖)、Wilms瘤/PMC11911335,TISCH2无对应数据集,继续搁置)。

HNSC_GSE139324与PMC12575298(候选池同义重复第3例)原文使用的数据集完全相同——可以用本研究
自己的z-score方法在同一份数据上复核该模型是否真的落在APOE+TAM对应的髓系/巨噬细胞类别,
作为候选池同义重复判断的额外经验验证(该模型无论结果如何都不计入H1a/H1b,因候选池本身
已被判定同义重复)。
"""
import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sp
from collections import Counter

OUT_DIR = "analysis_output/data/"
SDL = "/Users/yun/SignatureDL_raw_data/pan_cancer_phase2/"

TISCH2_JOBS = {
    "HNSC_GSE139324": {
        "dir": SDL + "HNSC_GSE139324/",
        "models": {
            "PMC12575298_NKTAM_HNSCC": ["ITGB7", "KDELR1", "ALG5", "ERP44", "BRI3", "TMBIM6",
                                          "BCAP31", "ATP6V0E1", "PLAU", "NDUFA4", "IGFLR1", "PSMD7",
                                          "PSMC1", "NDFIP1", "SNX6", "PSMB5", "TCEAL4", "ARL6IP1",
                                          "TRPV2", "P4HA1", "TPP1", "ATP2C1", "IDH1"],
            "PMC9945660_autophagy_lncRNA_HNSCC": ["AC010326.3", "AL160006.1", "AL122010.1",
                                                     "AC139530.1", "AC092747.4", "AL139287.1",
                                                     "MIR503HG", "AC009318.2", "LINC01711", "LINC02560"],
            "PMC8791745_hypoxia_lncRNA_HNSCC": ["AC116914.2", "AC144831.1", "AL357033.4",
                                                   "LINC00460", "LINC01980", "LINC02195", "MIAT",
                                                   "MSC-AS1", "MYOSLID"],
        },
    },
    "OSCC_GSE172577": {
        "dir": SDL + "OSCC_GSE172577/",
        "models": {
            "PMC12971589_CAF_OSCC": ["THBS1", "SLFN11", "TSPAN11", "ADAMTS15", "CCR7", "AQP1",
                                       "IL13RA2", "SERPINB7"],
            "PMC10599730_ERstress_OSCC": ["IBSP", "RDM1", "RBP4"],
        },
    },
    "ESCA_GSE160269": {
        "dir": SDL + "ESCA_GSE160269/",
        "models": {
            "PMC9271611_autophagy_lncRNA_ESCC": ["AC004690.2", "AC092159.3", "AC093627.4",
                                                    "AL078604.2", "BDNF-AS", "HAND2-AS1", "LINC00410",
                                                    "LINC00588", "PSMD6-AS2", "ZEB1-AS1", "LINC02586"],
            "PMC11603746_disulfidptosis_lncRNA_ESCC": ["NALT1", "AC118755.1", "LINC01770", "IPO5P1",
                                                          "AC104041.1", "AC138207.5", "AC092484.1",
                                                          "AC083799.1", "AC012467.2"],
        },
    },
    "UVM_GSE139829": {
        "dir": SDL + "UVM_GSE139829/",
        "models": {
            "PMC10287976_chromatin_UVM": ["RUVBL1", "SIRT3", "SMARCD3"],
            "PMC9709208_ICD_UVM": ["CASP8", "ENTPD1", "FOXP3", "IL6", "LY96"],
        },
    },
    "ALL_GSE132509": {
        "dir": SDL + "ALL_GSE132509/",
        "models": {
            "PMC12078231_ubiquitination_ALL": ["ATL2", "MKRN1", "FBXW8", "FBXO8", "DCAF16", "WSB1",
                                                  "CHFR", "MDM2", "SOCS2"],
        },
    },
    "OV_GSE154600": {
        "dir": SDL + "OV_GSE154600/",
        "models": {
            "PMC8464158_glycolysis_lncRNA_OV": ["AC133644.2", "CTD-2396E7.11", "CTD-3065J16.9",
                                                   "LINC00240", "TMEM254-AS1"],
        },
    },
    "PRAD_GSE172301": {
        "dir": SDL + "PRAD_GSE172301/",
        "models": {
            "PMC11550951_cuproptosis_lncRNA_PRAD": ["AC010896.1", "AC016394.2", "SNHG9"],
        },
    },
    "MB_GSE155446": {
        "dir": SDL + "MB_GSE155446/",
        "models": {
            "PMC12504172_nucleotide_MB": ["EIF4EBP1", "VARS1", "PRKDC", "KARS1", "SCARB1", "COL4A5",
                                            "SMN1", "FASN", "PLPBP", "PTPN1", "ATAD3A", "APTX",
                                            "EIF2AK3", "LYST", "PDE4D", "HDAC8", "PINK1"],
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
    ct_col = "Celltype (major-lineage)" if "Celltype (major-lineage)" in meta.columns else meta.columns[
        [c for c in meta.columns if "major" in c.lower() or "celltype" in c.lower()][0]
    ]
    common = expr_df.index.intersection(meta.index)
    expr_df = expr_df.loc[common]
    expr_df["celltype"] = meta.loc[common, ct_col].values
    mean_by_ct = expr_df.groupby("celltype")[present].mean().T
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)

    rows = []
    for model_name, raw_genes in models.items():
        mp = [resolved[g] for g in raw_genes if g in resolved]
        if not mp:
            print(f"\n=== {model_name} (0/{len(raw_genes)}) === 完全不可分析,跳过")
            continue
        sub = pd.DataFrame({"model": model_name, "gene": mp, "top_celltype": top_ct[mp].values,
                             "zscore": z.loc[mp].max(axis=1).round(2).values})
        rows.append(sub)
        print(f"\n=== {model_name} ({len(mp)}/{len(raw_genes)}) ===")
        print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
        print("分布:", dict(Counter(sub.top_celltype)))
        print("候选细胞类型(K):", len(mean_by_ct.columns), list(mean_by_ct.columns))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main():
    all_dfs = []
    for tag, cfg in TISCH2_JOBS.items():
        try:
            df = run_tisch2(tag, cfg)
        except Exception as e:
            print(f"\n########## {tag} 跳过(数据未就绪: {e}) ##########")
            continue
        if len(df):
            all_dfs.append(df)
    out = pd.concat(all_dfs, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_newcancertype_attribution.csv", index=False)
    print(f"\n写入 {OUT_DIR}phase2_newcancertype_attribution.csv, 总模型数: {out['model'].nunique()}, 总基因数: {len(out)}")
    print(f"模型清单: {sorted(out['model'].unique())}")


if __name__ == "__main__":
    main()
