"""
第5条任务:挖掘已复核的93个命中模型原文,寻找"作者自述的单细胞观察"作为独立、非同义重复的
佐证证据——即原作者自己对本模型基因做的单细胞/scRNA-seq分析中,报告了基因富集/表达于某个
细胞类型的观察(不论作者本人是否从"命名是否合理"这一角度解读),这类观察独立于本研究自己的
z-score归因流程(不同团队、可能不同数据集、不同分析代码),可以作为跨团队交叉印证使用,
参考前序血管生成课题3.3节"模型#8/#9"和3.9节"Zhou et al. 2025"的先例。

方法:对每个模型的PMC全文,先定位一个"单细胞/scRNA-seq"上下文窗口(段落级别),再检查该窗口
内是否同时出现该模型至少一个基因符号+一个细胞类型关键词。命中的候选逐一人工复核确认:
(a) 该单细胞分析确实是作者自己对本模型基因做的(不是引用他人研究、不是候选池构建阶段的
    marker基因鉴定——那属于候选池构建方法,不是"验证"意义上的观察);
(b) 报告的细胞类型是否是TME细胞(基质/免疫)而非命名对应的通路预期细胞类型,或反过来是否
    确认了命名对应的机制(两种情况都要如实记录,不是只找支持本研究假说的证据)。
"""
import re
from pathlib import Path
from collections import defaultdict

import pandas as pd

DATA_DIR = "analysis_output/data/"
PKG_DIR = Path("analysis_output/data/pmc_packages")

CELLTYPE_KEYWORDS = [
    "macrophage", "fibroblast", "CAF", "cancer-associated fibroblast", "endothelial",
    "T cell", "CD8", "CD4", "Treg", "B cell", "NK cell", "natural killer",
    "dendritic cell", "monocyte", "myeloid", "neutrophil", "mast cell",
    "pericyte", "smooth muscle", "epithelial cell", "malignant cell", "tumor cell",
    "stromal cell", "plasma cell", "TAM", "M2", "M1 polariz",
]
SC_CONTEXT_KEYWORDS = [
    "single-cell", "single cell", "scRNA-seq", "scRNA seq", "sc-RNA", "10x Genomics",
    "UMAP", "t-SNE", "tSNE", "cell cluster", "cell annotation",
]
CANDIDATE_POOL_EXCLUDE_HINTS = [
    "candidate gene", "candidate pool", "marker genes of", "module genes",
    "WGCNA", "FindAllMarkers", "differentially expressed genes between",
]

def load_model_genes():
    df = pd.read_csv(DATA_DIR + "phase2_pilot_all_combined.csv")
    model_genes = defaultdict(set)
    for _, row in df.iterrows():
        model_genes[row["model"]].add(row["gene"].upper())
    return model_genes

def get_fulltext(pmcid):
    pkg = PKG_DIR / pmcid
    if not pkg.exists():
        return None
    txt_files = sorted(pkg.glob(f"{pmcid}*.txt"))
    if not txt_files:
        return None
    return txt_files[0].read_text(encoding="utf-8", errors="ignore")

def find_sc_paragraphs(text):
    paras = re.split(r"\n\s*\n", text)
    hits = []
    for p in paras:
        if len(p) < 40:
            continue
        if any(kw.lower() in p.lower() for kw in SC_CONTEXT_KEYWORDS):
            hits.append(p)
    return hits

def scan_model(model, genes, text):
    paras = find_sc_paragraphs(text)
    findings = []
    gene_pattern = re.compile(
        r"\b(" + "|".join(re.escape(g) for g in genes if len(g) >= 2) + r")\b"
    )
    for p in paras:
        gene_matches = set(m.upper() for m in gene_pattern.findall(p))
        if not gene_matches:
            continue
        ct_matches = [kw for kw in CELLTYPE_KEYWORDS if kw.lower() in p.lower()]
        if not ct_matches:
            continue
        is_pool_construction = any(h.lower() in p.lower() for h in CANDIDATE_POOL_EXCLUDE_HINTS)
        findings.append({
            "model": model,
            "genes_mentioned": sorted(gene_matches),
            "celltypes_mentioned": ct_matches,
            "likely_pool_construction": is_pool_construction,
            "paragraph": p.strip()[:1200],
        })
    return findings

def main():
    model_genes = load_model_genes()
    all_findings = []
    n_scanned = 0
    n_no_fulltext = 0
    for model, genes in sorted(model_genes.items()):
        pmcid = model.split("_")[0]
        text = get_fulltext(pmcid)
        if text is None:
            n_no_fulltext += 1
            continue
        n_scanned += 1
        findings = scan_model(model, genes, text)
        all_findings.extend(findings)

    print(f"扫描模型数: {n_scanned}, 无全文可读: {n_no_fulltext}")
    print(f"候选段落命中数(基因+细胞类型关键词同段出现): {len(all_findings)}")
    non_pool = [f for f in all_findings if not f["likely_pool_construction"]]
    print(f"其中排除候选池构建上下文后剩余: {len(non_pool)}")

    out_rows = []
    for f in all_findings:
        out_rows.append({
            "model": f["model"],
            "genes_mentioned": ";".join(f["genes_mentioned"]),
            "celltypes_mentioned": ";".join(f["celltypes_mentioned"]),
            "likely_pool_construction": f["likely_pool_construction"],
            "paragraph_excerpt": f["paragraph"].replace("\n", " "),
        })
    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(DATA_DIR + "author_reported_scRNA_candidates.csv", index=False)
    print(f"写入 {DATA_DIR}author_reported_scRNA_candidates.csv,供逐条人工复核")

if __name__ == "__main__":
    main()
