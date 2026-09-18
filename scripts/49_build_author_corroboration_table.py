# -*- coding: utf-8 -*-
"""
把作者自述单细胞观察佐证(3.3/5节引用的`author_reported_scRNA_curated.json`)的7个
模型案例整理成表格,供审稿人直接查看逐基因比对结果,而不必解析JSON。

输出:
  analysis_output/data/author_reported_corroboration_table.csv/.md
"""
import json
from pathlib import Path

import pandas as pd


def main():
    d = json.load(open("analysis_output/data/author_reported_scRNA_curated.json", encoding="utf-8"))
    cases = d["cases"]
    df = pd.DataFrame(cases)
    df.to_csv("analysis_output/data/author_reported_corroboration_table.csv", index=False)

    lines = [
        "| # | 模型 | 主题 | 癌种 | 作者数据集/方法 | 作者发现 | 本研究数据集 | 本研究归因结果 | 一致性 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for i, c in enumerate(cases, 1):
        lines.append(
            f"| {i} | {c['model']} | {c['theme']} | {c['cancer']} | {c['author_dataset']} | "
            f"{c['author_finding']} | {c['our_dataset']} | {c['our_finding']} | {c['agreement']} |"
        )
    lines.append("")
    lines.append(f"**汇总**:{d['aggregate_summary']['n_models_with_usable_author_reported_observation']}个模型、"
                  f"{d['aggregate_summary']['n_gene_level_comparisons']}组基因级比对,"
                  f"{d['aggregate_summary']['n_exact_celltype_match']}组精确匹配、"
                  f"{d['aggregate_summary']['n_category_level_match_only']}组仅大类一致、"
                  f"{d['aggregate_summary']['n_mismatch']}组不一致。")
    lines.append("")
    lines.append(f"**不一致说明**:{d['aggregate_summary']['mismatch_note']}")
    Path("analysis_output/data/author_reported_corroboration_table.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )

    print(f"写入 author_reported_corroboration_table.csv/.md ({len(df)}行案例 + 汇总)")


if __name__ == "__main__":
    main()
