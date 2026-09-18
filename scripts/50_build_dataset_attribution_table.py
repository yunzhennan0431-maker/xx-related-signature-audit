# -*- coding: utf-8 -*-
"""
生成全部21个单细胞参考数据集的完整归因汇总表:谱系覆盖度分级(K值、覆盖大类数、
最小细胞类型样本量)与实际使用该数据集的模型数/显著模型数/H1a-H1b分布合并为一张表。

输出:
  analysis_output/data/dataset_attribution_summary_21.csv/.md
"""
from pathlib import Path

import pandas as pd

CANCER_NAMES = {
    "GSE81861": "结直肠癌", "GSE176078": "乳腺癌", "GSE151530": "肝细胞癌",
    "GSE167297": "胃癌", "SKCM_GSE120575": "黑色素瘤", "Glioma_GSE131928": "胶质瘤",
    "AML_GSE116256": "急性髓系白血病", "KIRC_GSE159115": "肾透明细胞癌",
    "OS_GSE162454": "骨肉瘤", "UCEC_GSE139555": "子宫内膜癌", "PAAD_GSE154778": "胰腺癌",
    "NSCLC_GSE131907": "非小细胞肺癌", "CESC_GSE168652": "宫颈癌", "BLCA_GSE130001": "膀胱癌",
    "HNSC_GSE139324": "头颈鳞癌", "OSCC_GSE172577": "口腔鳞癌", "UVM_GSE139829": "葡萄膜黑色素瘤",
    "ALL_GSE132509": "急性淋巴细胞白血病", "OV_GSE154600": "卵巢癌", "PRAD_GSE172301": "前列腺癌",
    "MB_GSE155446": "髓母细胞瘤",
}

H1A_CATS = {"Stromal_Vascular", "Immune_Hematopoietic"}
H1B_CATS = {"Malignant"}


def main():
    grading = pd.read_csv("analysis_output/data/dataset_quality_grading.csv")
    combined = pd.read_csv("analysis_output/data/phase3_meta_analysis_with_miss_pool_recovery.csv")

    agg_rows = []
    for ds, g in combined.groupby("dataset"):
        n_models = len(g)
        n_sig = int(g["significant_fdr05"].sum())
        n_h1a = int(((g["top_category"].isin(H1A_CATS)) & g["significant_fdr05"] & ~g["tautological_pool"]).sum())
        n_h1b = int(((g["top_category"].isin(H1B_CATS)) & g["significant_fdr05"] & ~g["tautological_pool"]).sum())
        agg_rows.append({
            "dataset": ds, "281模型中使用该数据集的模型数": n_models,
            "其中FDR显著数": n_sig, "H1a(TME集中,非同义重复)": n_h1a,
            "H1b(肿瘤自身集中,非同义重复)": n_h1b,
        })
    agg = pd.DataFrame(agg_rows)

    merged = grading.merge(agg, on="dataset", how="left").fillna(
        {"281模型中使用该数据集的模型数": 0, "其中FDR显著数": 0,
         "H1a(TME集中,非同义重复)": 0, "H1b(肿瘤自身集中,非同义重复)": 0}
    )
    for c in ["281模型中使用该数据集的模型数", "其中FDR显著数", "H1a(TME集中,非同义重复)", "H1b(肿瘤自身集中,非同义重复)"]:
        merged[c] = merged[c].astype(int)
    merged.insert(1, "癌种", merged["dataset"].map(CANCER_NAMES))
    merged = merged.sort_values("dataset").reset_index(drop=True)

    merged.to_csv("analysis_output/data/dataset_attribution_summary_21.csv", index=False)

    lines = [
        "| 数据集 | 癌种 | K(候选细胞类型数) | 谱系覆盖(0-4大类) | 使用该数据集的模型数(281中) | 显著模型数 | H1a | H1b |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for _, r in merged.iterrows():
        lines.append(
            f"| {r['dataset']} | {r['癌种']} | {r['K_celltypes']} | {r['lineage_coverage_0to4']} | "
            f"{r['281模型中使用该数据集的模型数']} | {r['其中FDR显著数']} | "
            f"{r['H1a(TME集中,非同义重复)']} | {r['H1b(肿瘤自身集中,非同义重复)']} |"
        )
    Path("analysis_output/data/dataset_attribution_summary_21.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )

    print(f"写入 dataset_attribution_summary_21.csv/.md ({len(merged)}个数据集)")
    print("覆盖度<4/4的数据集:", merged[merged["lineage_coverage_0to4"] < 4]["dataset"].tolist())


if __name__ == "__main__":
    main()
