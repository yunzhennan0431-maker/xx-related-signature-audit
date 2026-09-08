"""
系统性排查全部108个模型是否存在类似PMC10900655/PMC13114778的重复发表/误标癌种问题。
局限性第5节已承认这一排查此前只是"理论上可自动化但未执行"的待办,本脚本补上。

方法:两两比对模型的基因集合(标准化大小写、"."与"-"互换后取并集),计算Jaccard相似度。
基因集合完全重叠或高度重叠(且面板不太小,避免小面板偶然巧合)的模型对,标记为需要人工核查
是否也存在系数表重复。这不是充分条件(相同基因集不代表系数也相同),但足以作为初筛,
真正确认仍需人工回读原文系数表(参照PMC10900655/PMC13114778的核实方式)。
"""
import re
from itertools import combinations

import pandas as pd

DATA_DIR = "analysis_output/data/"


def normalize_gene(g):
    g = g.strip().upper()
    return g


def main():
    df = pd.read_csv(DATA_DIR + "phase2_pilot_all_combined.csv")
    df["gene_norm"] = df["gene"].apply(normalize_gene)
    model_genes = df.groupby("model")["gene_norm"].apply(lambda g: frozenset(g)).to_dict()
    models = sorted(model_genes.keys())
    print(f"共{len(models)}个模型,开始两两比对({len(models)*(len(models)-1)//2}对)")

    rows = []
    for m1, m2 in combinations(models, 2):
        g1, g2 = model_genes[m1], model_genes[m2]
        if len(g1) < 3 or len(g2) < 3:
            continue  # 面板过小(<3基因)时高重叠很可能是巧合,不纳入筛查
        inter = g1 & g2
        union = g1 | g2
        jaccard = len(inter) / len(union) if union else 0
        if jaccard >= 0.5:
            rows.append({
                "model1": m1, "model2": m2,
                "n1": len(g1), "n2": len(g2), "n_shared": len(inter),
                "jaccard": round(jaccard, 3),
                "shared_genes": ";".join(sorted(inter)),
            })

    cols = ["model1", "model2", "n1", "n2", "n_shared", "jaccard", "shared_genes"]
    if not rows:
        print("未发现Jaccard相似度>=0.5的模型对(面板>=3基因)")
        out = pd.DataFrame(columns=cols)
    else:
        out = pd.DataFrame(rows).sort_values("jaccard", ascending=False)
        print(f"\n发现{len(out)}对候选,按Jaccard相似度排序:")
        print(out[["model1", "model2", "n1", "n2", "n_shared", "jaccard"]].to_string(index=False))

    out.to_csv(DATA_DIR + "duplicate_model_scan_candidates.csv", index=False)
    print(f"\n写入 {DATA_DIR}duplicate_model_scan_candidates.csv"
          f"(即使为空也写入,留痕说明本次全量排查确实执行过、结果是零候选,而非未运行)")


if __name__ == "__main__":
    main()
