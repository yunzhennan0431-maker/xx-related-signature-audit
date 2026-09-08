# -*- coding: utf-8 -*-
"""
局限性第5节承认:33号脚本的重复发表排查只比对了基因集合(Jaccard相似度),没有比对原始
回归系数数值,因此无法排除"基因集合不同但仍存在部分系数篡改/重复使用"这类更隐蔽的问题。

本脚本补做一次真正的系数值级别比对——但只能覆盖仍保留原始系数的模型子集:batch1-3的
verified_extractions_batch{1,2,3}.json里,人工复核时把每个基因连同其回归系数一起记录了下来
(共27个模型);batch6-11改用了更精简的hits/misses结构,只记录基因数与数据集,没有保留系数值,
因此无法纳入这次比对。这不是"不完整的系数比对"的借口,而是如实说明:全量系数级排查需要重新
回读全部108个模型的原始文献表格,是一次接近重做部分抽取工作量级的任务,超出这次补做的范围。

方法:对27个模型两两比较共有基因的系数值(允许极小浮点误差,阈值1e-6),报告完全重合基因子集
里系数值精确匹配的比例。理由是:即使两个模型的基因集合不完全相同(例如一个是另一个的子集,或
只有部分重叠),如果重叠部分的系数值精确到小数点后多位仍然完全一致,那就和PMC10900655/
PMC13114778案例一样,是重复使用同一套系数表的强证据,不能仅凭基因集合不同就排除嫌疑。
"""
import json
import glob
from itertools import combinations

import pandas as pd

FILES = sorted(glob.glob("analysis_output/data/verified_extractions_batch1.json")) + \
        sorted(glob.glob("analysis_output/data/verified_extractions_batch2.json")) + \
        sorted(glob.glob("analysis_output/data/verified_extractions_batch3.json"))


def load_models():
    models = {}
    for fp in FILES:
        d = json.load(open(fp))
        for x in d.get("results", []):
            if not isinstance(x, dict) or x.get("verdict") != "hit":
                continue
            genes = x.get("genes", [])
            if len(genes) < 2 or not all("coef" in g for g in genes):
                continue
            pmcid = x.get("pmcid")
            models[pmcid] = {g["gene"].strip().upper(): g["coef"] for g in genes}
    return models


def main():
    models = load_models()
    print(f"可用于系数级比对的模型数: {len(models)}(仅batch1-3,保留了原始系数值)")
    ids = sorted(models.keys())
    rows = []
    for m1, m2 in combinations(ids, 2):
        g1, g2 = models[m1], models[m2]
        shared = set(g1) & set(g2)
        if len(shared) < 2:
            continue
        exact_match = sum(1 for g in shared if abs(g1[g] - g2[g]) < 1e-6)
        rows.append({
            "model1": m1, "model2": m2,
            "n1": len(g1), "n2": len(g2), "n_shared_genes": len(shared),
            "n_coef_exact_match": exact_match,
            "coef_match_rate": round(exact_match / len(shared), 3),
            "shared_genes": ";".join(sorted(shared)),
        })
    out = pd.DataFrame(rows, columns=["model1", "model2", "n1", "n2", "n_shared_genes",
                                       "n_coef_exact_match", "coef_match_rate", "shared_genes"])
    out = out.sort_values("coef_match_rate", ascending=False)
    out.to_csv("analysis_output/data/coefficient_level_duplicate_check_batch1to3.csv", index=False)

    print(f"\n共{len(ids)*(len(ids)-1)//2}对两两比较,其中{len(out)}对有>=2个共享基因")
    suspicious = out[out.coef_match_rate >= 0.8]
    print(f"\n系数值精确匹配率>=80%的可疑候选对: {len(suspicious)}")
    if len(suspicious):
        print(suspicious.to_string(index=False))
    else:
        print("(无)")
    print("\n共享基因数最多的几对(供人工复核参考,即使匹配率不高):")
    print(out.sort_values("n_shared_genes", ascending=False).head(5).to_string(index=False))


if __name__ == "__main__":
    main()
