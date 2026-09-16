"""
对"miss池系统性复查"发现的、疑似被原始复核错误判为miss的论文(见44号脚本同批产出的
miss_pool_recheck_results.csv/confident_hit_genes_v2.csv/ready_for_attribution.json),
补充跑z-score归因,评估把这些模型并入108模型集合后,论文headline结论(23.1%显著率、
16:6 H1a:H1b)会发生多大变化。

复用41号脚本已经写好的各数据集loader(log-normalize过的表达矩阵+细胞类型标注),
自己实现z-score归因(与项目一贯方法一致:按细胞类型求均值→跨细胞类型z-score→argmax),
不复用41号脚本的wilcoxon_attribution。
"""
import importlib.util
import json
import numpy as np
import pandas as pd

spec41 = importlib.util.spec_from_file_location("m41", "analysis_output/scripts/41_full_wilcoxon_robustness_check.py")
m41 = importlib.util.module_from_spec(spec41)
spec41.loader.exec_module(m41)

spec15 = importlib.util.spec_from_file_location("m15", "analysis_output/scripts/15_cross_theme_meta_analysis.py")
m15 = importlib.util.module_from_spec(spec15)
spec15.loader.exec_module(m15)

DATA_DIR = "analysis_output/data/"


def zscore_attribution(expr_df: pd.DataFrame, celltype: pd.Series, genes) -> dict:
    present = [g for g in genes if g in expr_df.columns]
    if not present:
        return {}
    mean_by_ct = {}
    for ct in celltype.unique():
        cells = celltype[celltype == ct].index
        mean_by_ct[ct] = expr_df.loc[cells, present].mean(axis=0)
    mean_df = pd.DataFrame(mean_by_ct)
    z = mean_df.sub(mean_df.mean(axis=1), axis=0).div(mean_df.std(axis=1) + 1e-9, axis=0)
    top_ct = z.idxmax(axis=1)
    return {g: top_ct[g] for g in present}


def main():
    ready = json.load(open("/tmp/ready_for_attribution_TRULY_final.json"))
    by_dataset = {}
    for r in ready:
        by_dataset.setdefault(r["dataset"], []).append(r)

    all_rows = []
    for dataset, papers in by_dataset.items():
        all_genes = sorted(set(g for p in papers for g in p["genes"].split(";") if g))
        print(f"[{dataset}] {len(papers)} papers, {len(all_genes)} unique genes to look up")
        try:
            if dataset == "GSE81861":
                expr, ct = m41.load_gse81861(all_genes)
            elif dataset == "GSE176078":
                expr, ct = m41.load_brca(all_genes)
            elif dataset == "GSE151530":
                expr, ct = m41.load_hcc(all_genes)
            elif dataset in m41.TISCH2_PATHS:
                expr, ct = m41.load_tisch2(m41.TISCH2_PATHS[dataset], all_genes)
            else:
                print(f"  !! no loader for {dataset}, skipping {len(papers)} papers")
                continue
        except Exception as e:
            print(f"  !! failed to load {dataset}: {e}")
            continue

        gene_to_ct = zscore_attribution(expr, ct, all_genes)
        print(f"  -> attributed {len(gene_to_ct)}/{len(all_genes)} genes")

        for p in papers:
            genes = [g for g in p["genes"].split(";") if g]
            n_attributed = 0
            for g in genes:
                if g in gene_to_ct:
                    all_rows.append({"model": p["pmcid"], "gene": g, "top_celltype": gene_to_ct[g]})
                    n_attributed += 1
            if n_attributed < len(genes):
                print(f"    {p['pmcid']}: only {n_attributed}/{len(genes)} genes found in reference atlas")

    df = pd.DataFrame(all_rows)
    df.to_csv(DATA_DIR + "miss_pool_recovery_attribution.csv", index=False)
    print(f"\nTotal: {df['model'].nunique()} models, {len(df)} gene-model rows written to miss_pool_recovery_attribution.csv")


if __name__ == "__main__":
    main()
