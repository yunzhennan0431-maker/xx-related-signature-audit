# -*- coding: utf-8 -*-
"""
PMC9092832(消化系统乏氧相关lncRNA签名,ESCA+STAD+LIHC+CRC四癌种联合建模)的
多图谱敏感性检验:原分析只借用STAD(胃癌)图谱作为代表性单细胞参考做归因,是
方法学上的权宜之计。本脚本改在本研究已有的STAD/HCC(肝癌)/CRC(结直肠癌)三个
图谱上分别独立重跑同一套z-score归因,检验归因结果是否对"借用哪个癌种图谱"敏感。
(第四个来源癌种ESCA因ESCA_GSE160269下载未成功,不在本次比较范围内。)
"""
import numpy as np
import pandas as pd
from scipy.io import mmread
import h5py
import scipy.sparse as sp
from collections import Counter

GENES = ["LUCAT1", "MIR4435-2HG", "LINC01711", "AP000695.2", "ADAMTS9-AS2", "AC087521.1"]

CRC_PAN = "/Users/yun/CRC_raw_data/pan_cancer/"
CRC81861 = "/Users/yun/CRC_raw_data/GSE81861/"


def run_tisch2_h5(tag, data_dir, genes):
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
    print(f"[{tag}] present: {len(present)}/{len(genes)}  missing: {[g for g in genes if g not in resolved]}")

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
    sub = pd.DataFrame({"atlas": tag, "gene": present, "top_celltype": top_ct[present].values,
                         "zscore": z.loc[present].max(axis=1).round(2).values})
    return sub


def run_hcc(genes):
    tag = "HCC_GSE151530"
    data_dir = CRC_PAN + "HCC_GSE151530/"
    gene_tbl = pd.read_csv(data_dir + "GSE151530_genes.tsv.gz", sep="\t", header=None, names=["ensembl", "symbol"])
    barcodes = pd.read_csv(data_dir + "GSE151530_barcodes.tsv.gz", header=None, names=["barcode"])["barcode"].tolist()
    info = pd.read_csv(data_dir + "info.txt", sep="\t").set_index("Cell")
    mat = mmread(data_dir + "GSE151530_matrix.mtx.gz").tocsr()
    name_to_idx = {}
    for i, n in enumerate(gene_tbl["symbol"]):
        name_to_idx.setdefault(n, []).append(i)
    present = [g for g in genes if g in name_to_idx]
    print(f"[{tag}] present: {len(present)}/{len(genes)}  missing: {[g for g in genes if g not in name_to_idx]}")

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
    sub = pd.DataFrame({"atlas": tag, "gene": present, "top_celltype": top_ct[present].values,
                         "zscore": z.loc[present].max(axis=1).round(2).values})
    return sub


def run_crc(genes):
    tag = "CRC_GSE81861"

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

    genes_present = [g for g in genes if g in combined.index]
    print(f"[{tag}] present: {len(genes_present)}/{len(genes)}  missing: {[g for g in genes if g not in combined.index]}")
    logexpr = np.log2(combined.loc[genes_present] + 1)
    result = {}
    for ct in meta.celltype.unique():
        cells = meta[meta.celltype == ct].cell.values
        result[ct] = logexpr[cells].mean(axis=1)
    mean_by_ct = pd.DataFrame(result)
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)
    sub = pd.DataFrame({"atlas": tag, "gene": genes_present, "top_celltype": top_ct[genes_present].values,
                         "zscore": z.loc[genes_present].max(axis=1).round(2).values})
    return sub


def main():
    stad = run_tisch2_h5("STAD_GSE167297", CRC_PAN + "STAD_GSE167297/", GENES)
    hcc = run_hcc(GENES)
    crc = run_crc(GENES)
    out = pd.concat([stad, hcc, crc], ignore_index=True)
    out.to_csv("analysis_output/data/pmc9092832_multi_atlas_sensitivity.csv", index=False)

    print("\n===== 三图谱归因结果对比(长表) =====")
    print(out.sort_values(["gene", "atlas"]).to_string(index=False))

    print("\n===== 按基因透视(每个基因在三个图谱下的众数细胞类型) =====")
    pivot = out.pivot(index="gene", columns="atlas", values="top_celltype")
    print(pivot.to_string())

    print("\n===== 三图谱之间归因细胞类型一致性 =====")
    for g in GENES:
        if g in pivot.index:
            vals = pivot.loc[g].dropna()
            print(f"{g}: {dict(vals)}  {'一致' if vals.nunique() == 1 else '不一致'}")
        else:
            print(f"{g}: 至少一个图谱缺失该基因")


if __name__ == "__main__":
    main()
