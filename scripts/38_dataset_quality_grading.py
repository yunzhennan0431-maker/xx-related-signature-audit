# -*- coding: utf-8 -*-
"""
局限性第5节承认:"当前统计框架仍未对'覆盖不完整'数据集做系统性的分级加权,只是做了排除法
的事后验证"。本脚本补上一版具体、可复现的轻量分级方案,而不是继续停留在"以后可以考虑"这句
空话上。

分级方案(三个独立维度,不合成单一分数——合成单一分数会掩盖问题的具体来源,三个维度分开报告
更利于读者自行判断):
1. 谱系覆盖完整度(lineage_coverage,0-4整数):候选细胞类型能覆盖本研究四大类
   (Malignant/Stromal_Vascular/Immune_Hematopoietic/Epithelial_mixed)中的几类。这是最直接
   反映"数据集是否可能系统性缺失某类细胞、导致归因被迫落入偏态候选池"的指标,3.5/5节已识别的
   两个问题数据集(NSCLC_GSE131907覆盖1类、BLCA_GSE130001覆盖1类)预期会在此维度得分垫底。
2. 候选细胞类型数(K,越大意味着归因的"随机基线"越低、真实集中越不容易被小K巧合稀释,但K过大
   也会引入更多小样本类型的z-score噪声,不是越大越好,只作为描述性指标)。
3. 最小候选类型细胞数(min_n,现有敏感性分析已用"细胞数<20排除"这一阈值处理估计噪声问题,这
   里报告原始最小值,让读者看到具体逼近或跌破这一阈值的程度,而不只是"是否被排除"这一二元结果)。

三个维度都直接从`15_cross_theme_meta_analysis.py`里已经维护的`DATASET_CELLTYPE_COUNTS`(21个
数据集的细胞类型-细胞数字典)和`CELLTYPE_CATEGORY`(细胞类型到四大类的映射)计算得到,不需要
重新访问原始h5/mtx文件。
"""
import importlib.util

import pandas as pd

spec = importlib.util.spec_from_file_location(
    "m15", "analysis_output/scripts/15_cross_theme_meta_analysis.py")
m15 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m15)

DATASET_CELLTYPE_COUNTS = m15.DATASET_CELLTYPE_COUNTS
CELLTYPE_CATEGORY = m15.CELLTYPE_CATEGORY

FOUR_CATEGORIES = ["Malignant", "Stromal_Vascular", "Immune_Hematopoietic", "Epithelial_mixed"]


def main():
    rows = []
    for ds, counts in DATASET_CELLTYPE_COUNTS.items():
        K = len(counts)
        total_n = sum(counts.values())
        min_n = min(counts.values())
        min_ct = min(counts, key=counts.get)
        cats_present = set()
        unmapped = []
        for ct in counts:
            cat = CELLTYPE_CATEGORY.get(ct)
            if cat is None:
                unmapped.append(ct)
            else:
                cats_present.add(cat)
        lineage_coverage = len(cats_present & set(FOUR_CATEGORIES))
        missing_cats = [c for c in FOUR_CATEGORIES if c not in cats_present]
        rows.append({
            "dataset": ds, "K_celltypes": K, "total_cells": total_n,
            "lineage_coverage_0to4": lineage_coverage,
            "missing_categories": ";".join(missing_cats) if missing_cats else "(无缺失)",
            "min_celltype_n": min_n, "min_celltype_name": min_ct,
            "below_20cell_sensitivity_threshold": min_n < 20,
            "unmapped_celltypes": ";".join(unmapped) if unmapped else "",
        })

    df = pd.DataFrame(rows).sort_values(["lineage_coverage_0to4", "K_celltypes"])
    df.to_csv("analysis_output/data/dataset_quality_grading.csv", index=False)

    print(f"共{len(df)}个数据集完成分级\n")
    print(df.to_string(index=False))

    print("\n=== 谱系覆盖度<4(存在系统性缺失大类)的数据集 ===")
    flagged = df[df.lineage_coverage_0to4 < 4]
    print(flagged[["dataset", "lineage_coverage_0to4", "missing_categories", "K_celltypes"]].to_string(index=False))

    print(f"\n=== 与既有3.5/5节'覆盖不完整'认定的一致性核查 ===")
    known_flagged = {"NSCLC_GSE131907", "BLCA_GSE130001"}
    detected_flagged = set(flagged["dataset"])
    print("既有认定:", known_flagged)
    print("本次分级同样标记为覆盖度<4:", detected_flagged)
    print("既有认定是否均被本次量化分级重新识别:", known_flagged.issubset(detected_flagged))
    print("本次分级新增识别、既有正文未单独讨论的:", detected_flagged - known_flagged)


if __name__ == "__main__":
    main()
