# -*- coding: utf-8 -*-
"""
batch11复核中发现的重要质量问题:PMC13114778(2026,乳腺癌"LRGS"乳酸化签名,批次10已计入
102模型集)的23基因系数表与PMC10900655(2024,结直肠癌"LRGS"乳酸化签名,更早发表于
J Transl Med)逐基因、逐系数完全一致(精确到小数点后15位),而两篇论文标题、癌种、期刊
均不同——几乎可以确定PMC13114778是对PMC10900655模型输出的重复使用/误标癌种(而非独立
推导的乳腺癌模型)。按本研究"每个模型必须是独立衍生的签名"这一基本假设,不能把两者当作
两个独立样本点。处理方式:PMC10900655发表更早,判定为原始/合法模型,保留计入;
PMC13114778判定为疑似重复/误标,从102模型集中移除,不再作为独立乳腺癌模型统计。

本脚本对PMC10900655的23个基因(与PMC13114778完全相同的基因列表),按其实际声称的癌种
(结直肠癌)用GSE81861重新跑z-score归因,替换掉此前错误的乳腺癌(GSE176078)归因结果。
"""
import numpy as np
import pandas as pd
from collections import Counter

CRC81861 = "/Users/yun/CRC_raw_data/GSE81861/"
OUT_DIR = "analysis_output/data/"

GENES = ["HSPA1B", "LGALS4", "HSPB1", "ARPC1B", "ACTG1", "LEPROTL1", "HMGN2", "SFPQ", "RBM3",
         "SLC2A3", "LY6E", "ARL4C", "UBE2I", "AP2M1", "TERF2IP", "H2AFY", "RAB5C", "LITAF",
         "MED10", "RBM17", "METTL9", "NR1H2", "PTTG1"]


def load(fp):
    df = pd.read_csv(fp, index_col=0)
    df.index = [str(i).split("_")[1] if len(str(i).split("_")) >= 3 else str(i) for i in df.index]
    return df


def main():
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

    genes_present = [g for g in GENES if g in combined.index]
    print("present:", len(genes_present), "/", len(GENES), "missing:", [g for g in GENES if g not in combined.index])
    logexpr = np.log2(combined.loc[genes_present] + 1)
    result = {}
    for ct in meta.celltype.unique():
        cells = meta[meta.celltype == ct].cell.values
        result[ct] = logexpr[cells].mean(axis=1)
    mean_by_ct = pd.DataFrame(result)
    z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)

    sub = pd.DataFrame({"model": "PMC10900655_lactylation_CRC", "gene": genes_present,
                         "top_celltype": top_ct[genes_present].values,
                         "zscore": z.loc[genes_present].max(axis=1).round(2).values})
    print(sub[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))
    print("分布:", dict(Counter(sub.top_celltype)))
    print("候选细胞类型(K):", len(mean_by_ct.columns), list(mean_by_ct.columns))
    sub.to_csv(OUT_DIR + "phase2_pmc10900655_crc_attribution.csv", index=False)
    print(f"\n写入 {OUT_DIR}phase2_pmc10900655_crc_attribution.csv")


if __name__ == "__main__":
    main()
