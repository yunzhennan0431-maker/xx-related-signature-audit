# -*- coding: utf-8 -*-
"""
把图2(显著集中比例11个checkpoint轨迹)与图3(跨主题异质性置换检验p值4个checkpoint
轨迹)背后的原始数字导出为表格,供审稿人核对精确数值(图上只标了小数点后1-3位)。
源数字与`34_trajectory_figures.py`中硬编码的绘图数据完全一致。

输出:
  analysis_output/data/significant_rate_trajectory.csv/.md   (11个checkpoint)
  analysis_output/data/heterogeneity_pvalue_trajectory.csv/.md (4个checkpoint)
"""
from pathlib import Path

import pandas as pd

# 与34_trajectory_figures.py中n_models_a/rates_a完全一致
N_MODELS_A = [19, 25, 27, 37, 48, 61, 76, 77, 93, 102, 108]
RATES_A = [47.4, 40.0, 33.3, 24.3, 20.8, 13.1, 19.7, 20.8, 20.4, 19.6, 23.1]

# 与34_trajectory_figures.py中n_models_b/pvals_b完全一致；
# 77模型取论文报告区间(0.19-0.21)的中点0.20作代表值，93/102/108为脚本实测精确值
N_MODELS_B = [77, 93, 102, 108]
PVALS_B = [0.20, 0.0787, 0.1547, 0.3220]
PVALS_B_NOTE = ["区间中点代表值(报告区间0.19-0.21)", "脚本实测精确值", "脚本实测精确值", "脚本实测精确值"]


def main():
    df_a = pd.DataFrame({
        "checkpoint序号": range(1, len(N_MODELS_A) + 1),
        "累计纳入模型数N": N_MODELS_A,
        "FDR小于0.05显著集中比例(%)": RATES_A,
    })
    df_a.to_csv("analysis_output/data/significant_rate_trajectory.csv", index=False)
    md_a = ["| checkpoint | N | 显著集中比例(%) |", "|---|---|---|"]
    for _, r in df_a.iterrows():
        md_a.append(f"| {r['checkpoint序号']} | {r['累计纳入模型数N']} | {r['FDR小于0.05显著集中比例(%)']:.1f} |")
    Path("analysis_output/data/significant_rate_trajectory.md").write_text("\n".join(md_a), encoding="utf-8")

    df_b = pd.DataFrame({
        "checkpoint序号": range(1, len(N_MODELS_B) + 1),
        "累计纳入模型数N": N_MODELS_B,
        "置换检验p值": PVALS_B,
        "说明": PVALS_B_NOTE,
    })
    df_b.to_csv("analysis_output/data/heterogeneity_pvalue_trajectory.csv", index=False)
    md_b = ["| checkpoint | N | p值 | 说明 |", "|---|---|---|---|"]
    for _, r in df_b.iterrows():
        md_b.append(f"| {r['checkpoint序号']} | {r['累计纳入模型数N']} | {r['置换检验p值']:.4f} | {r['说明']} |")
    Path("analysis_output/data/heterogeneity_pvalue_trajectory.md").write_text("\n".join(md_b), encoding="utf-8")

    print(f"写入 significant_rate_trajectory.csv/.md ({len(df_a)}行)")
    print(f"写入 heterogeneity_pvalue_trajectory.csv/.md ({len(df_b)}行)")


if __name__ == "__main__":
    main()
