# -*- coding: utf-8 -*-
"""
38号脚本(数据集质量分级)用"候选细胞类型覆盖本研究四大类中的几类"这一量化指标重新核查
21个数据集后发现,除此前已知的NSCLC_GSE131907、BLCA_GSE130001外,SKCM_GSE120575与
HNSC_GSE139324同样只覆盖1/4大类(候选类型全部是免疫细胞亚型,不含恶性/基质/上皮细胞)——
这意味着依赖这4个数据集的模型(11+7+2+5=25个),其归因结果在恶性细胞方向(H1b)上是
结构性不可能出现的,不管真实生物学信号如何,都只能落在TME细胞(H1a)一侧。

本脚本把15号脚本"敏感性分析2"的排除范围从{NSCLC_GSE131907, BLCA_GSE130001}扩大到
{NSCLC_GSE131907, BLCA_GSE130001, SKCM_GSE120575, HNSC_GSE139324},重新计算显著比例
与H1a/H1b分布,检验主结论对这一更大范围排除是否依然稳健。
"""
import importlib.util

import pandas as pd

spec = importlib.util.spec_from_file_location(
    "m15", "analysis_output/scripts/15_cross_theme_meta_analysis.py")
m15 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m15)

DATA_DIR = "analysis_output/data/"


def main():
    df = pd.read_csv(m15.COMBINED_PATH)
    main_res = m15.analyze(df, exclude_rare=False)
    main_sig = main_res[main_res["significant_fdr05"]]
    main_sig_genuine = main_sig[~main_sig["tautological_pool"]]
    main_h1a = main_sig_genuine[main_sig_genuine["top_category"].isin(["Stromal_Vascular", "Immune_Hematopoietic"])]
    main_h1b = main_sig_genuine[main_sig_genuine["top_category"] == "Malignant"]

    EXPANDED_INCOMPLETE_DATASETS = {"NSCLC_GSE131907", "BLCA_GSE130001", "SKCM_GSE120575", "HNSC_GSE139324"}
    excluded_models = {m for m, d in m15.MODEL_TO_DATASET.items() if d in EXPANDED_INCOMPLETE_DATASETS}
    df_cov = df[~df["model"].isin(excluded_models)]
    cov_res = m15.analyze(df_cov, exclude_rare=False)
    cov_res.to_csv(DATA_DIR + "phase3_meta_analysis_expanded_coverage_sensitivity.csv", index=False)

    print(f"排除{len(excluded_models)}个模型(来自{sorted(EXPANDED_INCOMPLETE_DATASETS)}),"
          f"剩余{len(cov_res)}个模型  (对比:原有敏感性分析2只排除16个模型,剩余92个)")
    print(f"显著模型数(FDR<0.05): {cov_res['significant_fdr05'].sum()}/{len(cov_res)} "
          f"({100*cov_res['significant_fdr05'].sum()/len(cov_res):.1f}%)  "
          f"(主分析为 {main_res['significant_fdr05'].sum()}/{len(main_res)} "
          f"({100*main_res['significant_fdr05'].sum()/len(main_res):.1f}%))")

    cov_sig = cov_res[cov_res["significant_fdr05"]]
    cov_sig_genuine = cov_sig[~cov_sig["tautological_pool"]]
    cov_h1a = cov_sig_genuine[cov_sig_genuine["top_category"].isin(["Stromal_Vascular", "Immune_Hematopoietic"])]
    cov_h1b = cov_sig_genuine[cov_sig_genuine["top_category"] == "Malignant"]
    print(f"剔除同义重复后,H1a={len(cov_h1a)}, H1b={len(cov_h1b)}  "
          f"(主分析为 H1a={len(main_h1a)}, H1b={len(main_h1b)})")

    print("\n被排除模型清单(按数据集分组):")
    for ds in sorted(EXPANDED_INCOMPLETE_DATASETS):
        ms = sorted(m for m, d in m15.MODEL_TO_DATASET.items() if d == ds)
        print(f"  {ds} ({len(ms)}个): {ms}")

    print(f"\n写入 {DATA_DIR}phase3_meta_analysis_expanded_coverage_sensitivity.csv")


if __name__ == "__main__":
    main()
