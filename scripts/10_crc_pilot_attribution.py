"""
阶段2试点:用前序课题("血管生成签名横向对比研究")已验证过的z-score单细胞归因方法
(参考实现: Nutstore/Identification on angiogenesis/analysis_output/scripts/analyze_ge81861.py),
对本课题阶段1抽取到的2个结直肠癌(CRC)signature模型做基因-细胞类型归因,验证自动化流水线
能否在新课题里跑通、结果是否合理。

方法与前序课题完全一致,不是"新方法",这一步只是换了输入基因列表、复用同一套验证过的逻辑,
作为阶段2更大规模自动化归因(SCimilarity等)之前的最简对照锚点。

数据: /Users/yun/CRC_raw_data/GSE81861/ (与前序课题共用,不重复下载)
"""
import pandas as pd
import numpy as np

DATA_DIR = "/Users/yun/CRC_raw_data/GSE81861/"
OUT_DIR = "analysis_output/data/"

MODELS = {
    "PMC13190656_ferroptosis_CRC": [
        "ALG3", "CST1", "RPS10", "TRPV4", "SPRR1A", "FGFR4", "TPM1", "KDM1A",
    ],
    "PMC11656000_pyroptosis_CRC": [
        "IL20RB", "MID2", "IFITM10", "LAMP5", "CALB2", "BANK1", "ANGPTL4", "CCL22",
        "EGFL7", "UPK3B", "SPTBN5", "TMPRSS11E", "LINGO1", "FENDRR", "TRAF1", "RNF207",
        "ROBO3", "TMEM88", "GRP", "SYNGR3", "HOXC11", "CHGB", "HEYL", "P2RX5", "TOX2",
    ],
    # 以下2个模型的原文癌种是"colon cancer"/"colon adenocarcinoma",和GSE81861的结直肠癌
    # (colorectal, 含结肠+直肠)高度重合,复用同一份数据,不额外找结肠癌专属图谱。
    "PMC11384320_m6A_ferroptosis_COAD": ["HSD17B11", "VEGFA", "CXCL2", "ASNS", "FABP4", "GPX2"],
    "PMC11297001_oxidative_stress_COAD": [
        "RYR2", "GSR", "GSTM1", "HSPA1A", "ACADL", "STK25", "ALOX12", "MAPK12",
        "SERPINA1", "CYP19A1", "NOL3", "NGF", "GDF15", "NTRK2", "DAPK1", "UCN", "HBA2",
    ],
}


def load(fp):
    df = pd.read_csv(fp, index_col=0)
    df.index = [str(i).split("_")[1] if len(str(i).split("_")) >= 3 else str(i) for i in df.index]
    return df


def main():
    nm = load(DATA_DIR + "GSE81861_CRC_NM_all_cells_FPKM.csv")
    tm = load(DATA_DIR + "GSE81861_CRC_tumor_all_cells_FPKM.csv")
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

    print(f"细胞类型分布(GSE81861全量): \n{meta.celltype.value_counts()}\n")

    all_summaries = []
    for model_name, genes in MODELS.items():
        genes_present = [g for g in genes if g in combined.index]
        genes_missing = [g for g in genes if g not in combined.index]
        if genes_missing:
            print(f"[{model_name}] 未在GSE81861里找到的基因(可能是lncRNA/曾用名/该数据集未测到): {genes_missing}")

        logexpr = np.log2(combined.loc[genes_present] + 1)
        result = {}
        for ct in meta.celltype.unique():
            cells = meta[meta.celltype == ct].cell.values
            result[ct] = logexpr[cells].mean(axis=1)
        mean_by_ct = pd.DataFrame(result)

        z = mean_by_ct.sub(mean_by_ct.mean(axis=1), axis=0).div(mean_by_ct.std(axis=1) + 1e-9, axis=0)
        top_ct = z.idxmax(axis=1)
        summary = pd.DataFrame({
            "model": model_name, "gene": top_ct.index, "top_celltype": top_ct.values,
            "zscore": z.max(axis=1).round(2).values,
        })
        all_summaries.append(summary)

        print(f"\n=== {model_name} ({len(genes_present)}/{len(genes)}个基因在数据集中找到) ===")
        print(summary[["gene", "top_celltype", "zscore"]].sort_values("top_celltype").to_string(index=False))

        from collections import Counter
        ct_dist = Counter(top_ct.values)
        print(f"\n细胞类型归因分布: {dict(ct_dist)}")

    out = pd.concat(all_summaries, ignore_index=True)
    out.to_csv(OUT_DIR + "phase2_crc_pilot_attribution.csv", index=False)
    print(f"\n写入 {OUT_DIR}phase2_crc_pilot_attribution.csv")


if __name__ == "__main__":
    main()
