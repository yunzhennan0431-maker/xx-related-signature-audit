# -*- coding: utf-8 -*-
"""
生成参考文献章节的附表:全部93个纳入meta分析的独立签名模型,列出PMCID、标题、癌种、
主题、基因数,供`论文初稿_中文.md`参考文献章节引用(附表形式,不在正文参考文献列表中逐条列出)。
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "analysis_output/scripts")
from importlib import import_module
meta = import_module("15_cross_theme_meta_analysis")

MODEL_TO_DATASET = meta.MODEL_TO_DATASET

DATASET_TO_CANCER = {
    "GSE81861": "结直肠癌",
    "GSE176078": "乳腺癌",
    "GSE151530": "肝细胞癌",
    "GSE167297": "胃癌",
    "SKCM_GSE120575": "黑色素瘤",
    "Glioma_GSE131928": "胶质瘤",
    "AML_GSE116256": "急性髓系白血病",
    "KIRC_GSE159115": "肾透明细胞癌",
    "OS_GSE162454": "骨肉瘤",
    "UCEC_GSE139555": "子宫内膜癌",
    "PAAD_GSE154778": "胰腺癌",
    "NSCLC_GSE131907": "非小细胞肺癌",
    "CESC_GSE168652": "宫颈癌",
    "BLCA_GSE130001": "膀胱癌",
    "HNSC_GSE139324": "头颈鳞癌",
    "OSCC_GSE172577": "口腔鳞癌",
    "UVM_GSE139829": "葡萄膜黑色素瘤",
    "ALL_GSE132509": "急性淋巴细胞白血病",
    "OV_GSE154600": "卵巢癌",
    "PRAD_GSE172301": "前列腺癌",
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
    genes = pd.read_csv("analysis_output/data/phase2_pilot_all_combined.csv")
    gene_counts = genes.groupby("model")["gene"].nunique().to_dict()

    rows = []
    for model in sorted(MODEL_TO_DATASET.keys()):
        pmcid = model.split("_")[0]
        theme_cancer_part = model[len(pmcid) + 1:]
        dataset = MODEL_TO_DATASET[model]
        cancer = DATASET_TO_CANCER.get(dataset, dataset)
        title = get_title(pmcid)
        n_genes = gene_counts.get(model, "")
        rows.append({
            "PMCID": pmcid,
            "标题": title,
            "癌种(归因参考图谱)": cancer,
            "模型标识(主题_癌种缩写)": theme_cancer_part,
            "基因数": n_genes,
        })

    df = pd.DataFrame(rows)
    df.to_csv("analysis_output/data/reference_table_allmodels.csv", index=False)
    print(f"共{len(df)}个模型,写入 analysis_output/data/reference_table_allmodels.csv")

    # 生成Markdown表格供直接粘贴进论文
    lines = ["| # | PMCID | 标题 | 癌种 | 主题标识 | 基因数 |", "|---|---|---|---|---|---|"]
    for i, row in enumerate(rows, 1):
        lines.append(
            f"| {i} | {row['PMCID']} | {row['标题']} | {row['癌种(归因参考图谱)']} | "
            f"{row['模型标识(主题_癌种缩写)']} | {row['基因数']} |"
        )
    Path("analysis_output/data/reference_table_allmodels.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    print("写入 analysis_output/data/reference_table_allmodels.md")

    missing = df[df["标题"] == "(标题未找到)"]
    if len(missing):
        print(f"警告: {len(missing)}个模型标题未找到: {missing['PMCID'].tolist()}")


if __name__ == "__main__":
    main()
