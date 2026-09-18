# -*- coding: utf-8 -*-
"""
生成合并后281个模型(原108 + miss池复查新增173)的完整清单,供审稿人/读者核对每个
模型的PMCID、标题、癌种、基因数、K值、归因结果、显著性、是否同义重复、来源批次。

输出:
  analysis_output/data/model_list_281.csv
  analysis_output/data/model_list_281.md (供直接粘贴的Markdown表)
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "analysis_output/scripts")

DATASET_TO_CANCER = {
    "GSE81861": "结直肠癌", "GSE176078": "乳腺癌", "GSE151530": "肝细胞癌",
    "GSE167297": "胃癌", "SKCM_GSE120575": "黑色素瘤", "Glioma_GSE131928": "胶质瘤",
    "AML_GSE116256": "急性髓系白血病", "KIRC_GSE159115": "肾透明细胞癌",
    "OS_GSE162454": "骨肉瘤", "UCEC_GSE139555": "子宫内膜癌", "PAAD_GSE154778": "胰腺癌",
    "NSCLC_GSE131907": "非小细胞肺癌", "CESC_GSE168652": "宫颈癌", "BLCA_GSE130001": "膀胱癌",
    "HNSC_GSE139324": "头颈鳞癌", "OSCC_GSE172577": "口腔鳞癌", "UVM_GSE139829": "葡萄膜黑色素瘤",
    "ALL_GSE132509": "急性淋巴细胞白血病", "OV_GSE154600": "卵巢癌", "PRAD_GSE172301": "前列腺癌",
    "MB_GSE155446": "髓母细胞瘤",
}

PKG_DIR = Path("analysis_output/data/pmc_packages")


def get_title(pmcid):
    pkg = PKG_DIR / pmcid
    if not pkg.exists():
        return "(标题未找到)"
    xml_files = sorted(pkg.glob(f"{pmcid}*.xml"))
    if not xml_files:
        return "(标题未找到)"
    text = xml_files[0].read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"<article-title>(.*?)</article-title>", text, re.DOTALL)
    if not m:
        return "(标题未找到)"
    title = re.sub(r"<[^>]+>", "", m.group(1))
    title = re.sub(r"\s+", " ", title).strip()
    return title


def main():
    combined = pd.read_csv("analysis_output/data/phase3_meta_analysis_with_miss_pool_recovery.csv")
    orig_main = pd.read_csv("analysis_output/data/phase3_meta_analysis_main.csv")
    orig_model_ids = set(orig_main["model"])  # exact full model_id, not just PMCID
    assert len(orig_model_ids) == 108, f"expected 108 original model_ids, got {len(orig_model_ids)}"

    rows = []
    for _, r in combined.iterrows():
        model = r["model"]
        pmcid = model.split("_")[0]
        source = "original_108" if model in orig_model_ids else "miss_pool_recovery_173"
        cancer = DATASET_TO_CANCER.get(r["dataset"], r["dataset"])
        title = get_title(pmcid)
        rows.append({
            "PMCID": pmcid,
            "model_id": model,
            "标题": title,
            "癌种(归因参考图谱)": cancer,
            "参考数据集": r["dataset"],
            "基因数": r["n_genes"],
            "K_细胞类型数": r["K_celltypes"],
            "众数细胞类型": r["top_celltype"],
            "众数细胞类型大类": r["top_category"],
            "二项检验p值": r["binom_pvalue"],
            "FDR校正q值": round(r["fdr_qvalue"], 6) if pd.notna(r["fdr_qvalue"]) else "",
            "FDR<0.05显著": r["significant_fdr05"],
            "基因数<6检验功效不足": r["underpowered_n_lt6"],
            "候选池同义重复": r["tautological_pool"],
            "来源批次": source,
        })

    df = pd.DataFrame(rows).sort_values(["来源批次", "PMCID"]).reset_index(drop=True)
    out_csv = Path("analysis_output/data/model_list_281.csv")
    df.to_csv(out_csv, index=False)
    print(f"共{len(df)}个模型(原108: {(df['来源批次']=='original_108').sum()}, "
          f"miss池新增173: {(df['来源批次']=='miss_pool_recovery_173').sum()}),写入 {out_csv}")

    lines = ["| # | PMCID | 标题 | 癌种 | 数据集 | 基因数 | K | 众数细胞类型 | 大类 | q值 | 显著 | 同义重复 | 来源 |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, row in enumerate(df.itertuples(index=False), 1):
        lines.append(
            f"| {i} | {row[0]} | {row[2]} | {row[3]} | {row[4]} | {row[5]} | {row[6]} | "
            f"{row[7]} | {row[8]} | {row[10]} | {row[11]} | {row[13]} | {row[14]} |"
        )
    out_md = Path("analysis_output/data/model_list_281.md")
    out_md.write_text("\n".join(lines), encoding="utf-8")
    print(f"写入 {out_md}")

    missing = df[df["标题"] == "(标题未找到)"]
    if len(missing):
        print(f"警告: {len(missing)}个模型标题未找到: {missing['PMCID'].tolist()}")


if __name__ == "__main__":
    main()
