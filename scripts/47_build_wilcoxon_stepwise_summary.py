# -*- coding: utf-8 -*-
"""
汇总独立算法(Wilcoxon)验证三步的逐步结果,供审稿人一次性查看每一步的样本量、
指标、与z-score的比对结果,而不必分别打开5个原始CSV。

输出:
  analysis_output/data/wilcoxon_stepwise_summary.csv/.md  (三步汇总表)
  analysis_output/data/wilcoxon_full_model_comparison_combined.csv  (108模型×z/Wilcoxon/下采样Wilcoxon 宽表)
"""
from pathlib import Path

import pandas as pd


def main():
    # ---- Step 1: 119组验证锚点比对 ----
    step1 = pd.read_csv("analysis_output/data/independent_wilcoxon_vs_zscore_broadcategory.csv")
    step1_n = len(step1)
    step1_exact = int(step1["match_vs_gold"].sum())
    step1_cat = int(step1["cat_match"].sum())

    # ---- Step 2: 全量108模型重新归因 ----
    step2 = pd.read_csv("analysis_output/data/wilcoxon_vs_zscore_full_model_comparison.csv")
    step2_n = len(step2)
    z_sig = int(step2["significant_fdr05_z"].sum())
    w_sig = int(step2["significant_fdr05_w"].sum())
    sig_agree_2 = int(step2["sig_agree"].sum())

    # ---- Step 3: 样本量下采样修正 ----
    step3 = pd.read_csv("analysis_output/data/wilcoxon_downsampled_vs_zscore_full_model_comparison.csv")
    step3_n = len(step3)
    ds_sig = int(step3["significant_fdr05_ds"].sum())
    sig_agree_3 = int(step3["sig_agree"].sum())

    summary_rows = [
        {
            "步骤": "第一步:验证锚点比对",
            "对象": "119组基因×数据集(与前序血管生成课题金标准)",
            "样本量": step1_n,
            "关键指标": "精确细胞类型一致率 / 三大类一致率",
            "结果": f"{step1_exact}/{step1_n} ({step1_exact/step1_n:.1%}) / {step1_cat}/{step1_n} ({step1_cat/step1_n:.1%})",
            "对应脚本": "40_independent_wilcoxon_validation.py",
            "对应数据文件": "independent_wilcoxon_vs_zscore_broadcategory.csv",
        },
        {
            "步骤": "第二步:全量108模型重新归因",
            "对象": "108个模型、21个数据集",
            "样本量": step2_n,
            "关键指标": "显著模型数(z-score vs Wilcoxon) / 模型级别显著性一致率",
            "结果": f"z={z_sig}({z_sig/step2_n:.1%}), w={w_sig}({w_sig/step2_n:.1%}) / {sig_agree_2}/{step2_n} ({sig_agree_2/step2_n:.1%})",
            "对应脚本": "41_full_wilcoxon_robustness_check.py",
            "对应数据文件": "wilcoxon_vs_zscore_full_model_comparison.csv",
        },
        {
            "步骤": "第三步:细胞类型样本量下采样修正",
            "对象": "108个模型(占优细胞类型下采样至组内中位数后重新归因)",
            "样本量": step3_n,
            "关键指标": "下采样后显著模型数 / 模型级别显著性一致率(相对z-score)",
            "结果": f"{ds_sig}/{step3_n} ({ds_sig/step3_n:.1%}) / {sig_agree_3}/{step3_n} ({sig_agree_3/step3_n:.1%})",
            "对应脚本": "42_full_wilcoxon_downsampled_check.py",
            "对应数据文件": "wilcoxon_downsampled_vs_zscore_full_model_comparison.csv",
        },
    ]
    summary = pd.DataFrame(summary_rows)
    summary.to_csv("analysis_output/data/wilcoxon_stepwise_summary.csv", index=False)

    lines = ["| 步骤 | 对象 | 样本量 | 关键指标 | 结果 | 对应脚本 |", "|---|---|---|---|---|---|"]
    for row in summary_rows:
        lines.append(
            f"| {row['步骤']} | {row['对象']} | {row['样本量']} | {row['关键指标']} | "
            f"{row['结果']} | `{row['对应脚本']}` |"
        )
    Path("analysis_output/data/wilcoxon_stepwise_summary.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )

    # ---- 合并的108模型宽表(z-score / Wilcoxon / 下采样Wilcoxon 并列) ----
    merged = step2.merge(
        step3[["model", "significant_fdr05_ds", "top_category_ds", "sig_agree"]]
        .rename(columns={"sig_agree": "sig_agree_ds_vs_z"}),
        on="model", how="outer",
    )
    merged = merged.rename(columns={"sig_agree": "sig_agree_w_vs_z"})
    merged.to_csv("analysis_output/data/wilcoxon_full_model_comparison_combined.csv", index=False)

    print(f"写入 wilcoxon_stepwise_summary.csv/.md ({len(summary)}行)")
    print(f"写入 wilcoxon_full_model_comparison_combined.csv ({len(merged)}行, {merged.shape[1]}列)")


if __name__ == "__main__":
    main()
