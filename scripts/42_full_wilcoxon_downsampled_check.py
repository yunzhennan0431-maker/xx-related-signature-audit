"""
直接检验41号脚本发现的"Wilcoxon系统性偏向样本量大的细胞类型(尤其是恶性细胞)"这一假设:
对全部21个数据集,把每个细胞类型的细胞数下采样到不超过该数据集细胞类型数的中位数
(高于中位数的类型被下采样到中位数,低于中位数的类型不动),消除掉组间样本量不均衡这个
混淆因素,再重新跑一遍Wilcoxon归因+全套统计,看显著率和H1a:H1b比值是否回到z-score结果
附近。这是对41号脚本发现的"重大偏离"最直接的机制检验,而不是停留在一个未经证实的推测上。

复用40/41号脚本的数据加载与Wilcoxon归因逻辑,以及15号脚本的统计分析逻辑不做修改。
"""
import importlib.util
import numpy as np
import pandas as pd

spec41 = importlib.util.spec_from_file_location("m41", "analysis_output/scripts/41_full_wilcoxon_robustness_check.py")
m41 = importlib.util.module_from_spec(spec41)
spec41.loader.exec_module(m41)

spec15 = importlib.util.spec_from_file_location("m15", "analysis_output/scripts/15_cross_theme_meta_analysis.py")
m15 = importlib.util.module_from_spec(spec15)
spec15.loader.exec_module(m15)

DATA_DIR = "analysis_output/data/"
RNG = np.random.RandomState(0)


def downsample_to_median(expr: pd.DataFrame, celltype: pd.Series):
    counts = celltype.value_counts()
    cap = int(np.median(counts.values))
    keep_idx = []
    for ct, n in counts.items():
        idx = celltype.index[celltype == ct]
        if n > cap:
            idx = RNG.choice(idx, size=cap, replace=False)
        keep_idx.extend(idx)
    keep_idx = pd.Index(keep_idx)
    return expr.loc[keep_idx], celltype.loc[keep_idx], cap, dict(counts)


def main():
    zscore_df = pd.read_csv(DATA_DIR + "phase2_pilot_all_combined.csv")
    zscore_df["dataset"] = zscore_df["model"].map(m15.MODEL_TO_DATASET)
    zscore_df = zscore_df.dropna(subset=["dataset"])

    all_rows = []
    for dataset, sub in zscore_df.groupby("dataset"):
        genes_needed = sorted(sub["gene"].unique())
        n_models = sub["model"].nunique()
        try:
            if dataset == "GSE81861":
                expr, ct = m41.load_gse81861(genes_needed)
            elif dataset == "GSE176078":
                expr, ct = m41.load_brca(genes_needed)
            elif dataset == "GSE151530":
                expr, ct = m41.load_hcc(genes_needed)
            elif dataset in m41.TISCH2_PATHS:
                expr, ct = m41.load_tisch2(m41.TISCH2_PATHS[dataset], genes_needed)
            else:
                print(f"  !! no loader for {dataset}, skipping")
                continue
        except Exception as e:
            print(f"  !! failed to load {dataset}: {e}")
            continue

        expr_ds, ct_ds, cap, before_counts = downsample_to_median(expr, ct)
        after_counts = dict(ct_ds.value_counts())
        print(f"[{dataset}] {n_models} models, {len(genes_needed)} genes, cap={cap}")
        print(f"   before: {before_counts}")
        print(f"   after:  {after_counts}")

        gene_to_ct = m41.wilcoxon_attribution(expr_ds, ct_ds, genes_needed)
        for _, row in sub.iterrows():
            ct_call = gene_to_ct.get(row["gene"])
            if ct_call is not None:
                all_rows.append({"model": row["model"], "gene": row["gene"], "top_celltype": ct_call})

    ds_df = pd.DataFrame(all_rows)
    ds_df.to_csv(DATA_DIR + "phase2_pilot_all_combined_wilcoxon_downsampled.csv", index=False)
    print(f"\nDownsampled Wilcoxon attribution covers {ds_df['model'].nunique()}/108 models, "
          f"{len(ds_df)}/{len(zscore_df)} gene-model rows")

    z_res = m15.analyze(zscore_df, exclude_rare=False)
    z_sig = z_res[z_res.significant_fdr05]
    z_sig_genuine = z_sig[~z_sig.tautological_pool]
    z_h1a = z_sig_genuine[z_sig_genuine.top_category.isin(["Stromal_Vascular", "Immune_Hematopoietic"])]
    z_h1b = z_sig_genuine[z_sig_genuine.top_category == "Malignant"]

    w_full = pd.read_csv(DATA_DIR + "phase2_pilot_all_combined_wilcoxon.csv")
    w_full_res = m15.analyze(w_full, exclude_rare=False)
    w_full_sig = w_full_res[w_full_res.significant_fdr05]
    w_full_genuine = w_full_sig[~w_full_sig.tautological_pool]
    w_full_h1a = w_full_genuine[w_full_genuine.top_category.isin(["Stromal_Vascular", "Immune_Hematopoietic"])]
    w_full_h1b = w_full_genuine[w_full_genuine.top_category == "Malignant"]

    ds_res = m15.analyze(ds_df, exclude_rare=False)
    ds_sig = ds_res[ds_res.significant_fdr05]
    ds_sig_genuine = ds_sig[~ds_sig.tautological_pool]
    ds_h1a = ds_sig_genuine[ds_sig_genuine.top_category.isin(["Stromal_Vascular", "Immune_Hematopoietic"])]
    ds_h1b = ds_sig_genuine[ds_sig_genuine.top_category == "Malignant"]
    ds_res.to_csv(DATA_DIR + "phase3_meta_analysis_wilcoxon_downsampled.csv", index=False)

    print("\n" + "=" * 70)
    print("三方对比")
    print("=" * 70)
    print(f"z-score(论文现有结论):              {len(z_sig)}/{len(z_res)}显著({len(z_sig)/len(z_res)*100:.1f}%), H1a={len(z_h1a)}, H1b={len(z_h1b)}")
    print(f"Wilcoxon(全量细胞,未下采样):         {len(w_full_sig)}/{len(w_full_res)}显著({len(w_full_sig)/len(w_full_res)*100:.1f}%), H1a={len(w_full_h1a)}, H1b={len(w_full_h1b)}")
    print(f"Wilcoxon(下采样至组间样本量均衡):     {len(ds_sig)}/{len(ds_res)}显著({len(ds_sig)/len(ds_res)*100:.1f}%), H1a={len(ds_h1a)}, H1b={len(ds_h1b)}")

    merged = z_res[["model", "significant_fdr05", "top_category"]].merge(
        ds_res[["model", "significant_fdr05", "top_category"]], on="model", suffixes=("_z", "_ds"), how="outer")
    merged["sig_agree"] = merged["significant_fdr05_z"] == merged["significant_fdr05_ds"]
    print(f"\n下采样Wilcoxon vs z-score, 模型级别'是否显著'一致: {merged['sig_agree'].sum()}/{len(merged)} "
          f"({merged['sig_agree'].sum()/len(merged)*100:.1f}%)")
    merged.to_csv(DATA_DIR + "wilcoxon_downsampled_vs_zscore_full_model_comparison.csv", index=False)


if __name__ == "__main__":
    main()
